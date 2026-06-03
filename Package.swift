// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "ReactorSim",
    platforms: [
        .macOS("26.0")          // macOS Tahoe — native Liquid Glass APIs
    ],
    targets: [
        .executableTarget(
            name: "ReactorSim",
            path: "Sources/ReactorSim"
        )
    ]
)
