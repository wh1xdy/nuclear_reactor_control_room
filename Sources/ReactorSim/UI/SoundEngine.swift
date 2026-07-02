// SoundEngine.swift — procedural control-room audio (no assets).
//
// Three voices, all synthesized in one AVAudioSourceNode render callback:
//  • Turbine hum — 50 Hz shaft fundamental + harmonics + a little broadband
//    hiss, amplitude and pitch tracking rpm (a coastdown audibly winds down).
//  • Annunciator chime — two short 880 Hz pulses when a NEW alarm comes in.
//  • Breaker clunk — a 60 ms filtered-noise thud on switchyard operation.
//
// Parameter writes come from the main thread (tiny torn-read risk is harmless
// for audio params); the render thread only reads. Master volume is modest —
// atmosphere, not a sound effects showcase.

import AVFoundation

final class SoundEngine: @unchecked Sendable {
    private let engine = AVAudioEngine()
    private var node: AVAudioSourceNode?
    private var started = false

    // Control parameters (main thread writes, render thread reads).
    private var humLevel: Double = 0        // 0…1 (rpm fraction)
    private var humPitch: Double = 1        // rpm fraction → harmonic pitch scale
    private var chimeAt:  Double = -1       // sample time of chime start (−1 = idle)
    private var clunkAt:  Double = -1

    var enabled: Bool = true { didSet { enabled ? start() : stop() } }

    private var t: Double = 0               // render-thread sample clock [s]
    private var noiseState: UInt64 = 0x9E3779B97F4A7C15

    func start() {
        guard !started, enabled else { return }
        let format = engine.outputNode.outputFormat(forBus: 0)
        let sr = format.sampleRate > 0 ? format.sampleRate : 44_100
        let src = AVAudioSourceNode { [weak self] _, _, frameCount, audioBufferList -> OSStatus in
            guard let self else { return noErr }
            let abl = UnsafeMutableAudioBufferListPointer(audioBufferList)
            let dt = 1.0 / sr
            for frame in 0..<Int(frameCount) {
                var s = 0.0
                // ── Turbine hum ──
                let lvl = self.humLevel
                if lvl > 0.01 {
                    let f0 = 50.0 * max(0.05, self.humPitch)
                    let w = 2 * Double.pi * self.t
                    s += 0.030 * lvl * sin(w * f0)
                    s += 0.018 * lvl * sin(w * f0 * 2)
                    s += 0.010 * lvl * sin(w * f0 * 3.02)
                    // broadband hiss (xorshift noise, cheap + deterministic)
                    self.noiseState ^= self.noiseState << 13
                    self.noiseState ^= self.noiseState >> 7
                    self.noiseState ^= self.noiseState << 17
                    let n = Double(Int64(bitPattern: self.noiseState)) / Double(Int64.max)
                    s += 0.008 * lvl * n
                }
                // ── Annunciator chime: two 880 Hz pulses ──
                if self.chimeAt >= 0 {
                    let a = self.t - self.chimeAt
                    if a >= 0, a < 0.9 {
                        let pulse = a < 0.35 ? a : (a >= 0.45 ? a - 0.45 : -1)
                        if pulse >= 0, pulse < 0.35 {
                            s += 0.10 * exp(-pulse * 9) * sin(2 * .pi * 880 * self.t)
                        }
                    } else if a >= 0.9 { self.chimeAt = -1 }
                }
                // ── Breaker clunk: short low noise burst ──
                if self.clunkAt >= 0 {
                    let a = self.t - self.clunkAt
                    if a >= 0, a < 0.09 {
                        self.noiseState ^= self.noiseState << 13
                        self.noiseState ^= self.noiseState >> 7
                        self.noiseState ^= self.noiseState << 17
                        let n = Double(Int64(bitPattern: self.noiseState)) / Double(Int64.max)
                        s += 0.16 * exp(-a * 45) * n * sin(2 * .pi * 90 * self.t)
                    } else if a >= 0.09 { self.clunkAt = -1 }
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

    func chime() { if started { chimeAt = t } }
    func clunk() { if started { clunkAt = t } }
}
