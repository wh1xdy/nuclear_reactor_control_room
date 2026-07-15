#!/bin/zsh
# make-app.sh — wrap the SwiftPM release binary in a proper ReactorSim.app.
# Usage: scripts/make-app.sh [path/to/icon-1024.png]
# Unique reverse-DNS bundle id (com.wh1xdy.ReactorSim): the generic "ReactorSim"
# id collided with a macOS Screen Time app limit keyed on that name, which
# popped a "Time Limit" dialog on launch. A distinct id Screen Time has never
# seen launches clean. The app defaults (skin/type) fall back to guided/PWR on
# first run — which is exactly the demo view — so losing the old defaults
# domain costs nothing.
set -euo pipefail
cd "$(dirname "$0")/.."

ICON_SRC="${1:-assets/icon-1024.png}"
APP=dist/ReactorSim.app

swift build -c release

rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"

cp .build/release/ReactorSim "$APP/Contents/MacOS/ReactorSim"
# SwiftPM resource bundle (WAV callouts, calibration.json) — the generated
# Bundle.module accessor searches Bundle.main.resourceURL.
cp -R .build/release/ReactorSim_ReactorSim.bundle "$APP/Contents/Resources/"

if [[ -f "$ICON_SRC" ]]; then
  ICONSET=$(mktemp -d)/AppIcon.iconset
  mkdir -p "$ICONSET"
  for s in 16 32 128 256 512; do
    sips -z $s $s             "$ICON_SRC" --out "$ICONSET/icon_${s}x${s}.png"      >/dev/null
    sips -z $((s*2)) $((s*2)) "$ICON_SRC" --out "$ICONSET/icon_${s}x${s}@2x.png"   >/dev/null
  done
  iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/AppIcon.icns"
fi

cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
	<key>CFBundleName</key>            <string>ReactorSim</string>
	<key>CFBundleDisplayName</key>     <string>ReactorSim</string>
	<key>CFBundleIdentifier</key>      <string>com.wh1xdy.ReactorSim</string>
	<key>CFBundleExecutable</key>      <string>ReactorSim</string>
	<key>CFBundleIconFile</key>        <string>AppIcon</string>
	<key>CFBundlePackageType</key>     <string>APPL</string>
	<key>CFBundleShortVersionString</key> <string>1.0</string>
	<key>CFBundleVersion</key>         <string>1</string>
	<key>LSMinimumSystemVersion</key>  <string>26.0</string>
	<key>NSHighResolutionCapable</key> <true/>
	<key>NSHumanReadableCopyright</key> <string>© 2026 wh1xdy · training simulator, not for real-plant use</string>
</dict>
</plist>
PLIST

codesign --force --deep -s - "$APP"

# Install to ~/Applications so Spotlight/Dock launch the fresh build (not a
# stale copy). Overwrite in place; re-register with Launch Services.
mkdir -p "$HOME/Applications"
rm -rf "$HOME/Applications/ReactorSim.app"
cp -R "$APP" "$HOME/Applications/ReactorSim.app"
/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister \
  -f "$HOME/Applications/ReactorSim.app" 2>/dev/null || true

echo "Built $APP and installed to ~/Applications/ReactorSim.app"
