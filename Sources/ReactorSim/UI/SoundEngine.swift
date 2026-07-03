// SoundEngine.swift — procedural control-room audio (no assets).
//
// Voices, all synthesized in one AVAudioSourceNode render callback:
//  • Turbine hum — 50 Hz shaft fundamental + harmonics + broadband hiss,
//    amplitude and pitch tracking rpm (a coastdown audibly winds down).
//  • Annunciator chime — two short 880 Hz pulses on an ordinary new alarm.
//  • Reactor-trip horn — URGENT alternating two-tone (620/470 Hz) for ~4 s;
//    a scram must not sound polite.
//  • Breaker operation — industrial double impact: deep ~55 Hz body with
//    falling pitch, a metallic snap transient, and the spring-mechanism
//    rebound ~70 ms later. Not a beep.
//  plus a SPOKEN annunciator (AVSpeechSynthesizer, calm female voice —
//  think No Man's Sky ship computer) for major events and, optionally,
//  every alarm.
//
// Pausing the simulation silences everything (output gate + speech pause).
// Parameter writes come from the main thread (tiny torn-read risk is
// harmless for audio params); the render thread only reads.

import AVFoundation

final class SoundEngine: @unchecked Sendable {
    private let engine = AVAudioEngine()
    private var node: AVAudioSourceNode?
    private var started = false
    private let speech = AVSpeechSynthesizer()
    private lazy var voice: AVSpeechSynthesisVoice? = {
        // Calm female en voice: prefer Samantha, else any enhanced en female,
        // else the en-US default (female on stock macOS).
        let all = AVSpeechSynthesisVoice.speechVoices()
        return all.first { $0.name.contains("Samantha") }
            ?? all.first { $0.language.hasPrefix("en") && $0.gender == .female }
            ?? AVSpeechSynthesisVoice(language: "en-US")
    }()

    // Control parameters (main thread writes, render thread reads).
    private var humLevel: Double = 0
    private var humPitch: Double = 1
    private var chimeAt:  Double = -1
    private var clunkAt:  Double = -1
    private var hornAt:   Double = -1
    private var gateTarget: Double = 1      // 0 while the sim is paused
    private var gate:       Double = 1      // ramped toward the target (~10 ms, no step click)

    /// Off under XCTest — `swift test` must not play klaxons and speak.
    var enabled: Bool = NSClassFromString("XCTestCase") == nil {
        didSet { enabled ? start() : stop() }
    }

    private var t: Double = 0               // render-thread sample clock [s]
    private var noiseState: UInt64 = 0x9E3779B97F4A7C15

    private func noise() -> Double {
        noiseState ^= noiseState << 13
        noiseState ^= noiseState >> 7
        noiseState ^= noiseState << 17
        return Double(Int64(bitPattern: noiseState)) / Double(Int64.max)
    }

    func start() {
        guard !started, enabled else { return }
        let format = engine.outputNode.outputFormat(forBus: 0)
        let sr = format.sampleRate > 0 ? format.sampleRate : 44_100
        let src = AVAudioSourceNode { [weak self] _, _, frameCount, audioBufferList -> OSStatus in
            guard let self else { return noErr }
            let abl = UnsafeMutableAudioBufferListPointer(audioBufferList)
            let dt = 1.0 / sr
            for frame in 0..<Int(frameCount) {
                // Soft pause gate: ramp toward the target (~10 ms) so pausing
                // mid-horn doesn't click; once silent, FREEZE the sample clock
                // so one-shot timelines (horn/chime/clunk) resume where paused.
                self.gate += (self.gateTarget - self.gate) * 0.002
                if self.gate < 0.005, self.gateTarget == 0 {
                    for buf in abl { buf.mData?.assumingMemoryBound(to: Float.self)[frame] = 0 }
                    continue
                }
                var s = 0.0
                let g = self.gate
                // ── Turbine hum ──
                let lvl = self.humLevel * g
                if lvl > 0.01 {
                    let f0 = 50.0 * max(0.05, self.humPitch)
                    let w = 2 * Double.pi * self.t
                    s += 0.030 * lvl * sin(w * f0)
                    s += 0.018 * lvl * sin(w * f0 * 2)
                    s += 0.010 * lvl * sin(w * f0 * 3.02)
                    s += 0.008 * lvl * self.noise()
                }
                // ── Annunciator chime: two 880 Hz pulses ──
                if self.hornAt < 0, self.chimeAt >= 0 {   // horn overrides chime
                    let a = self.t - self.chimeAt
                    if a >= 0, a < 0.9 {
                        let pulse = a < 0.35 ? a : (a >= 0.45 ? a - 0.45 : -1)
                        if pulse >= 0, pulse < 0.35 {
                            // pulse-relative phase: every pulse starts at a zero
                            // crossing — no onset click under the full envelope.
                            s += 0.10 * g * exp(-pulse * 9) * sin(2 * .pi * 880 * pulse)
                        }
                    } else if a >= 0.9 { self.chimeAt = -1 }
                }
                // ── Reactor-trip horn: alternating two-tone, ~4 s, urgent ──
                if self.hornAt >= 0 {
                    let a = self.t - self.hornAt
                    if a >= 0, a < 4.2 {
                        let f = (Int(a / 0.30) % 2 == 0) ? 620.0 : 470.0
                        // Horn-relative phase: 620 and 470 Hz both complete an
                        // integer cycle count per 0.30 s segment, so phase is 0
                        // at every tone switch — no pop 13× per trip.
                        let base = sin(2 * .pi * f * a)
                        let h3   = 0.33 * sin(2 * .pi * f * 3 * a)
                        let env  = a < 0.02 ? a / 0.02 : (a > 3.9 ? max(0, (4.2 - a) / 0.3) : 1)
                        s += 0.16 * g * env * (base + h3)
                    } else if a >= 4.2 { self.hornAt = -1 }
                }
                // ── Breaker: industrial double impact ──
                if self.clunkAt >= 0 {
                    let a = self.t - self.clunkAt
                    if a >= 0, a < 0.40 {
                        // Impact 1: deep body, pitch falling 75→45 Hz.
                        if a < 0.16 {
                            let f = 75.0 - 190.0 * a
                            s += 0.34 * g * exp(-a * 26) * sin(2 * .pi * f * a)
                        }
                        // Metallic snap transient (first 9 ms, bright noise).
                        if a < 0.009 {
                            s += 0.22 * g * exp(-a * 320) * self.noise()
                        }
                        // Impact 2: spring-mechanism rebound at +70 ms, softer/lower.
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
        engine.mainMixerNode.outputVolume = 0.8
        node = src
        do { try engine.start(); started = true } catch { started = false }
    }

    func stop() {
        guard started else { return }
        speech.stopSpeaking(at: .immediate)
        engine.stop()
        if let n = node { engine.detach(n) }
        node = nil
        started = false
    }

    /// Called each physics step: rpm fraction drives hum level + pitch.
    func update(rpmFraction: Double) {
        humLevel = max(0, min(1, rpmFraction))
        humPitch = max(0, min(1.1, rpmFraction))
    }

    /// Sim pause: everything goes quiet (soft ramp, then the sample clock
    /// freezes so a mid-horn pause resumes exactly where it left off).
    func setPaused(_ paused: Bool) {
        gateTarget = paused ? 0 : 1
        if paused { speech.pauseSpeaking(at: .word) } else { speech.continueSpeaking() }
    }

    // The horn owns the annunciator channel: a chime raised during it is
    // dropped (not deferred into an orphaned fragment), and a horn cancels
    // any chime already sounding.
    func chime() { if started, hornAt < 0 { chimeAt = t } }
    func clunk() { if started { clunkAt = t } }
    func horn()  { if started { hornAt = t; chimeAt = -1 } }

    /// Spoken annunciator — calm, deliberate, female. Alarm text arrives in
    /// board shorthand; speak() lowercases it so the synthesizer doesn't spell
    /// out capitals, and expands the common abbreviations. `priority` flushes
    /// the queue first: a safety callout must never wait behind alarm chatter.
    func speak(_ text: String, priority: Bool = false) {
        guard started, enabled else { return }
        if priority, speech.isSpeaking { speech.stopSpeaking(at: .word) }
        var msg = text.lowercased()
        for (abbr, full) in [("rcs", "R C S"), ("pzr", "pressurizer"), ("t-avg", "T average"),
                             ("porv", "P O R V"), ("eccs", "E C C S"), ("sg", "steam generator")] {
            msg = msg.replacingOccurrences(of: abbr, with: full)
        }
        let u = AVSpeechUtterance(string: msg)
        u.voice = voice
        u.rate = 0.45
        u.pitchMultiplier = 0.98
        u.volume = 0.75
        u.preUtteranceDelay = 0.1
        speech.speak(u)
    }
}
