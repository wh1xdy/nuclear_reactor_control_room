// SoundEngine.swift — control-room audio.
//
// Three layers:
//  • Procedural bed (one AVAudioSourceNode): the 50 Hz machine-hall BURR
//    (dense harmonic stack with slight detune beating + slow amplitude
//    wobble — a transformer hum, not a thin sine), annunciator chime, and
//    the urgent reactor-trip horn. All phase-continuous; the pause gate
//    ramps ~10 ms then freezes the sample clock so one-shots resume intact.
//  • Sampled breaker: the real recording (sliced into close / open halves,
//    bundled as resources) — falls back to the procedural double-thud if
//    the resource is missing.
//  • Helmet voice — the No-Man's-Sky-style suit annunciator: the system
//    female voice rendered OFFLINE via AVSpeechSynthesizer.write, then
//    played through a processing graph:
//      player → time-pitch (−60 ¢, slightly synthetic)
//             → EQ (high-pass 220 Hz, +5 dB resonance at 1.5 kHz,
//                   low-pass ~7.8 kHz — the "inside the helmet" band)
//             → distortion (.speechRadioTower, low wet — digital crispness
//                   with a faint ring-mod character)
//             → reverb (.smallRoom, ~10 % — a 20–50 ms speaker-by-the-ear
//                   space, never a hall)
//             → dynamics compressor (every word even and clear)
//    Own FIFO queue; priority callouts flush it. Pause parks the player.
//
// swift test stays silent (engine disabled under XCTest).

import AVFoundation

final class SoundEngine: @unchecked Sendable {
    private let engine = AVAudioEngine()
    private var node: AVAudioSourceNode?
    private var started = false

    // ── Voice graph ─────────────────────────────────────────────────────────
    private let synth = AVSpeechSynthesizer()
    private let voicePlayer = AVAudioPlayerNode()
    private let voicePitch  = AVAudioUnitTimePitch()
    private let voiceEQ     = AVAudioUnitEQ(numberOfBands: 3)
    private let voiceDist   = AVAudioUnitDistortion()
    private let voiceVerb   = AVAudioUnitReverb()
    private let voiceComp: AVAudioUnitEffect = {
        var desc = AudioComponentDescription()
        desc.componentType = kAudioUnitType_Effect
        desc.componentSubType = kAudioUnitSubType_DynamicsProcessor
        desc.componentManufacturer = kAudioUnitManufacturer_Apple
        return AVAudioUnitEffect(audioComponentDescription: desc)
    }()
    private var voiceFormat: AVAudioFormat?
    private var voiceQueue: [String] = []
    private var voiceBusy = false
    private lazy var voice: AVSpeechSynthesisVoice? = {
        let all = AVSpeechSynthesisVoice.speechVoices()
        return all.first { $0.name.contains("Samantha") }
            ?? all.first { $0.language.hasPrefix("en") && $0.gender == .female }
            ?? AVSpeechSynthesisVoice(language: "en-US")
    }()

    // ── Breaker samples ─────────────────────────────────────────────────────
    private var breakerClose: AVAudioPlayer?
    private var breakerOpen:  AVAudioPlayer?

    // ── Pre-rendered PA callouts (Samantha through the sox PA chain; the
    //    future AI-clone voice replaces these WAVs file-for-file) ────────────
    private var calloutPlayer: AVAudioPlayer?
    private var calloutQueue: [URL] = []
    private var calloutGen = 0            // bumps on flush → stale completions no-op

    // ── Procedural-bed parameters (main writes, render reads) ───────────────
    private var humLevel: Double = 0
    private var humPitch: Double = 1
    private var dieselMode: Int = 0        // 0 off · 1 cranking · 2 running
    private var chimeAt:  Double = -1
    private var clunkAt:  Double = -1        // procedural fallback only
    private var hornAt:   Double = -1
    private var gateTarget: Double = 1
    private var gate:       Double = 1

    /// Off under XCTest — `swift test` must not play klaxons and speak.
    var enabled: Bool = NSClassFromString("XCTestCase") == nil {
        didSet { enabled ? start() : stop() }
    }

    private var t: Double = 0
    private var noiseState: UInt64 = 0x9E3779B97F4A7C15
    private func noise() -> Double {
        noiseState ^= noiseState << 13
        noiseState ^= noiseState >> 7
        noiseState ^= noiseState << 17
        return Double(Int64(bitPattern: noiseState)) / Double(Int64.max)
    }

    // MARK: — Lifecycle

    func start() {
        guard !started, enabled else { return }
        let format = engine.outputNode.outputFormat(forBus: 0)
        let sr = format.sampleRate > 0 ? format.sampleRate : 44_100

        let src = AVAudioSourceNode { [weak self] _, _, frameCount, audioBufferList -> OSStatus in
            guard let self else { return noErr }
            let abl = UnsafeMutableAudioBufferListPointer(audioBufferList)
            let dt = 1.0 / sr
            for frame in 0..<Int(frameCount) {
                self.gate += (self.gateTarget - self.gate) * 0.002
                if self.gate < 0.005, self.gateTarget == 0 {
                    for buf in abl { buf.mData?.assumingMemoryBound(to: Float.self)[frame] = 0 }
                    continue        // clock frozen: one-shots resume where paused
                }
                var s = 0.0
                let g = self.gate
                // ── 50 Hz machine-hall burr ──
                let lvl = self.humLevel * g
                if lvl > 0.01 {
                    let f0 = 50.0 * max(0.05, self.humPitch)
                    let w = 2 * Double.pi * self.t
                    // Slow wobble (0.7 Hz, ±18 %) — the "alive" quality.
                    let am = 1.0 + 0.18 * sin(w * 0.7)
                    // Dense harmonic stack, mid harmonics carrying the energy
                    // (laptop speakers barely reproduce 50 Hz itself), plus a
                    // detuned second fundamental for the slow beat.
                    var hum = 0.014 * sin(w * f0)
                    hum += 0.030 * sin(w * f0 * 2)
                    hum += 0.026 * sin(w * f0 * 3.004)
                    hum += 0.014 * sin(w * f0 * 4)
                    hum += 0.008 * sin(w * f0 * 5.01)
                    hum += 0.012 * sin(w * (f0 + 0.35))       // beat partner
                    hum += 0.007 * self.noise()                // breath
                    s += lvl * am * hum
                }
                // ── Emergency diesel generators ──
                if self.dieselMode == 1 {
                    // Cranking: slow starter chug — noise bursts at ~6.5 Hz.
                    let chug = max(0, sin(2 * .pi * 6.5 * self.t))
                    s += 0.11 * g * chug * chug * chug * self.noise()
                    s += 0.02 * g * sin(2 * .pi * 18 * self.t)
                } else if self.dieselMode == 2 {
                    // Running: 24.5 Hz firing fundamental + harmonics + breath.
                    let w = 2 * Double.pi * self.t
                    s += 0.030 * g * sin(w * 24.5)
                    s += 0.020 * g * sin(w * 49.2)
                    s += 0.010 * g * sin(w * 73.6)
                    s += 0.010 * g * self.noise()
                }
                // ── Annunciator chime (horn owns the channel) ──
                if self.hornAt < 0, self.chimeAt >= 0 {
                    let a = self.t - self.chimeAt
                    if a >= 0, a < 0.9 {
                        let pulse = a < 0.35 ? a : (a >= 0.45 ? a - 0.45 : -1)
                        if pulse >= 0, pulse < 0.35 {
                            s += 0.10 * g * exp(-pulse * 9) * sin(2 * .pi * 880 * pulse)
                        }
                    } else if a >= 0.9 { self.chimeAt = -1 }
                }
                // ── Reactor-trip horn ──
                if self.hornAt >= 0 {
                    let a = self.t - self.hornAt
                    if a >= 0, a < 4.2 {
                        let f = (Int(a / 0.30) % 2 == 0) ? 620.0 : 470.0
                        let base = sin(2 * .pi * f * a)        // integer cycles per segment → no pops
                        let h3   = 0.33 * sin(2 * .pi * f * 3 * a)
                        let env  = a < 0.02 ? a / 0.02 : (a > 3.9 ? max(0, (4.2 - a) / 0.3) : 1)
                        s += 0.16 * g * env * (base + h3)
                    } else if a >= 4.2 { self.hornAt = -1 }
                }
                // ── Breaker (procedural fallback when no sample bundled) ──
                if self.clunkAt >= 0 {
                    let a = self.t - self.clunkAt
                    if a >= 0, a < 0.40 {
                        if a < 0.16 {
                            let f = 75.0 - 190.0 * a
                            s += 0.34 * g * exp(-a * 26) * sin(2 * .pi * f * a)
                        }
                        if a < 0.009 { s += 0.22 * g * exp(-a * 320) * self.noise() }
                        let b = a - 0.070
                        if b >= 0, b < 0.20 {
                            let f = 58.0 - 90.0 * b
                            s += 0.20 * g * exp(-b * 22) * sin(2 * .pi * f * b)
                            if b < 0.006 { s += 0.10 * g * exp(-b * 350) * self.noise() }
                        }
                    } else if a >= 0.40 { self.clunkAt = -1 }
                }
                let v = Float(max(-1, min(1, s)))
                for buf in abl {
                    buf.mData?.assumingMemoryBound(to: Float.self)[frame] = v
                }
                self.t += dt
            }
            return noErr
        }
        engine.attach(src)
        engine.connect(src, to: engine.mainMixerNode, format: format)

        // Helmet-voice graph. Wired ONCE here, in STEREO at the engine's own
        // rate: the DynamicsProcessor AU rejects mono (kAudioUnitErr −10868),
        // and connecting at the synthesizer's 22 kHz buffer format throws an
        // uncaught NSException — the rendered speech is CONVERTED into this
        // fixed graph format instead (see playRendered).
        engine.attach(voicePlayer)
        engine.attach(voicePitch)
        engine.attach(voiceEQ)
        engine.attach(voiceDist)
        engine.attach(voiceVerb)
        engine.attach(voiceComp)
        voicePitch.pitch = -60                              // slightly synthetic
        let hp = voiceEQ.bands[0]
        hp.filterType = .highPass; hp.frequency = 220; hp.bypass = false
        let res = voiceEQ.bands[1]
        res.filterType = .parametric; res.frequency = 1500; res.bandwidth = 1.0
        res.gain = 5; res.bypass = false                    // helmet resonance
        let lp = voiceEQ.bands[2]
        lp.filterType = .lowPass; lp.frequency = 7800; lp.bypass = false
        voiceDist.loadFactoryPreset(.speechRadioTower)      // crisp + faint ring-mod
        voiceDist.wetDryMix = 14
        voiceVerb.loadFactoryPreset(.smallRoom)             // 20–50 ms, by the ear
        voiceVerb.wetDryMix = 10
        if let vf = AVAudioFormat(standardFormatWithSampleRate: sr, channels: 2) {
            voiceFormat = vf
            engine.connect(voicePlayer, to: voicePitch, format: vf)
            engine.connect(voicePitch,  to: voiceEQ,    format: vf)
            engine.connect(voiceEQ,     to: voiceDist,  format: vf)
            engine.connect(voiceDist,   to: voiceVerb,  format: vf)
            engine.connect(voiceVerb,   to: voiceComp,  format: vf)
            engine.connect(voiceComp,   to: engine.mainMixerNode, format: vf)
        }

        engine.mainMixerNode.outputVolume = 0.8
        node = src
        do { try engine.start(); started = true } catch { started = false }

        // Breaker samples from the bundle (nil → procedural fallback).
        breakerClose = (try? AVAudioPlayer(contentsOf:
            Bundle.module.url(forResource: "breaker-close", withExtension: "wav") ?? URL(fileURLWithPath: "/nonexistent")))
        breakerOpen = (try? AVAudioPlayer(contentsOf:
            Bundle.module.url(forResource: "breaker-open", withExtension: "wav") ?? URL(fileURLWithPath: "/nonexistent")))
        breakerClose?.volume = 0.55
        breakerOpen?.volume = 0.55
        breakerClose?.prepareToPlay()
        breakerOpen?.prepareToPlay()
    }

    func stop() {
        guard started else { return }
        voiceQueue.removeAll(); voiceBusy = false
        voicePlayer.stop()
        engine.stop()
        if let n = node { engine.detach(n) }
        node = nil
        started = false
    }

    // MARK: — Controls

    func update(rpmFraction: Double) {
        humLevel = max(0, min(1, rpmFraction))
        humPitch = max(0, min(1.1, rpmFraction))
    }

    /// Emergency diesel state: 0 off, 1 cranking, 2 running.
    func setDiesel(_ mode: Int) { dieselMode = mode }

    func setPaused(_ paused: Bool) {
        gateTarget = paused ? 0 : 1
        if paused {
            if voicePlayer.isPlaying { voicePlayer.pause() }
            breakerClose?.pause(); breakerOpen?.pause()
            calloutPlayer?.pause()
        } else {
            if started, voiceBusy { voicePlayer.play() }
            if let c = calloutPlayer, c.currentTime > 0 { c.play() }
        }
    }

    func chime() { if started, hornAt < 0 { chimeAt = t } }
    func horn()  { if started { hornAt = t; chimeAt = -1 } }

    /// Breaker operation: play the real recording for the matching direction
    /// (close = energize, open = de-energize); procedural thud as fallback.
    func clunk(closing: Bool = true) {
        guard started else { return }
        if let p = closing ? breakerClose : breakerOpen {
            p.currentTime = 0
            p.play()
        } else {
            clunkAt = t
        }
    }

    // MARK: — PA callouts (pre-rendered) + helmet-voice fallback

    /// Announce by KEY: plays the bundled pre-rendered PA sample if present
    /// (exact sound guaranteed), else falls back to live helmet-chain speech.
    /// `priority` flushes both the sample queue and the live-voice queue.
    func announce(key: String, text: String, priority: Bool = false) {
        guard started, enabled else { return }
        if let url = Bundle.module.url(forResource: "callout-\(key)", withExtension: "wav") {
            if priority {
                calloutGen += 1
                calloutQueue.removeAll()
                calloutPlayer?.stop(); calloutPlayer = nil
                voiceQueue.removeAll(); voicePlayer.stop(); voiceBusy = false
            }
            calloutQueue.append(url)
            pumpCalloutQueue()
        } else {
            speak(text, priority: priority)
        }
    }

    private func pumpCalloutQueue() {
        guard started, calloutPlayer == nil, !calloutQueue.isEmpty else { return }
        let url = calloutQueue.removeFirst()
        guard let p = try? AVAudioPlayer(contentsOf: url) else { pumpCalloutQueue(); return }
        p.volume = 0.85
        calloutPlayer = p
        if gateTarget > 0 { p.play() }
        let gen = calloutGen
        DispatchQueue.main.asyncAfter(deadline: .now() + p.duration + 0.25) { [weak self] in
            guard let self, self.calloutGen == gen else { return }
            self.calloutPlayer = nil
            self.pumpCalloutQueue()
        }
    }

    // MARK: — Helmet voice (live synthesis — voice-every-alarm fallback path)

    /// Suit-annunciator speech: slow, even, female, processed through the
    /// helmet chain. `priority` flushes anything queued or playing first.
    func speak(_ text: String, priority: Bool = false) {
        guard started, enabled else { return }
        var msg = text.lowercased()
        for (abbr, full) in [("rcs", "R C S"), ("pzr", "pressurizer"), ("t-avg", "T average"),
                             ("porv", "P O R V"), ("eccs", "E C C S"), ("sg", "steam generator")] {
            msg = msg.replacingOccurrences(of: abbr, with: full)
        }
        if priority {
            voiceQueue.removeAll()
            voicePlayer.stop()
            voiceBusy = false
        }
        voiceQueue.append(msg)
        pumpVoiceQueue()
    }

    private func pumpVoiceQueue() {
        guard started, !voiceBusy, !voiceQueue.isEmpty else { return }
        voiceBusy = true
        let msg = voiceQueue.removeFirst()

        let u = AVSpeechUtterance(string: msg)
        u.voice = voice
        u.rate = 0.40                      // pretty slow, deliberate
        u.pitchMultiplier = 0.98           // near-monotone delivery
        u.volume = 1.0

        // Render offline, collect the PCM, then play through the helmet chain.
        var buffers: [AVAudioPCMBuffer] = []
        synth.write(u) { [weak self] buffer in
            guard let self else { return }
            if let pcm = buffer as? AVAudioPCMBuffer, pcm.frameLength > 0 {
                buffers.append(pcm)
            } else {
                // Zero-length buffer = utterance finished rendering.
                DispatchQueue.main.async { self.playRendered(buffers) }
            }
        }
    }

    private func playRendered(_ buffers: [AVAudioPCMBuffer]) {
        guard started, engine.isRunning, let vf = voiceFormat,
              let inFmt = buffers.first?.format else {
            voiceBusy = false
            pumpVoiceQueue()
            return
        }
        // Concatenate the synthesizer's buffers (22 kHz mono), then convert
        // into the FIXED graph format (hw-rate stereo) in one pass. Never
        // reconnect the chain per utterance — that path threw NSExceptions.
        let totalIn = buffers.reduce(AVAudioFrameCount(0)) { $0 + $1.frameLength }
        guard totalIn > 0,
              let big = AVAudioPCMBuffer(pcmFormat: inFmt, frameCapacity: totalIn),
              let conv = AVAudioConverter(from: inFmt, to: vf) else {
            voiceBusy = false
            pumpVoiceQueue()
            return
        }
        for b in buffers {
            let n = Int(b.frameLength)
            for ch in 0..<Int(inFmt.channelCount) {
                memcpy(big.floatChannelData![ch].advanced(by: Int(big.frameLength)),
                       b.floatChannelData![ch], n * MemoryLayout<Float>.size)
            }
            big.frameLength += b.frameLength
        }
        let outCap = AVAudioFrameCount(Double(totalIn) * (vf.sampleRate / inFmt.sampleRate)) + 4096
        guard let out = AVAudioPCMBuffer(pcmFormat: vf, frameCapacity: outCap) else {
            voiceBusy = false
            pumpVoiceQueue()
            return
        }
        var fed = false
        var err: NSError?
        let status = conv.convert(to: out, error: &err) { _, st in
            if fed { st.pointee = .endOfStream; return nil }
            fed = true; st.pointee = .haveData; return big
        }
        guard status != .error, out.frameLength > 0 else {
            voiceBusy = false
            pumpVoiceQueue()
            return
        }
        voicePlayer.scheduleBuffer(out) { [weak self] in
            DispatchQueue.main.async {
                self?.voiceBusy = false
                self?.pumpVoiceQueue()
            }
        }
        if gateTarget > 0 { voicePlayer.play() }
    }
}
