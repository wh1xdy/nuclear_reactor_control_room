# Nuclear Reactor Control Room Sim

En enkel kontrollrumssimulator i Python/Pygame med tre reaktormodeller:

- PWR
- BWR
- RBMK

## Kom igång

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
nrcr
```

Alternativt:

```bash
python run_control_room.py
```

## Styrning i simulatorn

- `1/2/3`: välj reaktortyp (PWR/BWR/RBMK)
- `W/S`: styrstavar
- `A/D`: flöde/pumphastighet
- `Q/E`: turbinventil
- `SHIFT`: snabb justering
- `SPACE`: SCRAM
- `R`: reset
- `ESC`: avsluta

## macOS universal `.dmg`

Repo:t har nu ett CI-flöde som bygger två `.app`-bundles (x86_64 + arm64), slår ihop dem till en universal `.app` och paketerar en DMG-artefakt.

Kvar för en full distributionskedja för slutanvändare:

1. **Kodsignering** av `.app` och `.dmg` (`codesign`, Developer ID).
2. **Notarisering** + stapling (`notarytool` + `stapler`).
3. **Releasepublicering** (t.ex. bifoga DMG till GitHub Release).
## Vad som saknas för en fungerande **universell** `.dmg` (macOS)

Just nu innehåller repot själva Python-koden, men ingen macOS-packetering. För en fungerande universal `.dmg` saknas i praktiken:

1. **Byggsteg för `.app`-bundle**
   - Lägg till t.ex. `pyinstaller` eller `py2app` i build-flödet.
   - Projektet har i nuläget bara Python-paketmetadata i `pyproject.toml`, inte app-bundle-konfiguration.

2. **Universal2-bygge (arm64 + x86_64)**
   - Bygg på ett sätt som ger universal2-binär (eller bygg två varianter och slå ihop med `lipo` där det behövs).
   - Säkerställ att Python-runtime och native beroenden (exempelvis `pygame`/SDL) också är universal2.

3. **DMG-skapande**
   - Lägg till steg/script för att skapa `.dmg` från `.app` (t.ex. `create-dmg`, `hdiutil` eller motsvarande CI-script).

4. **Kodsignering**
   - Signera `.app` och `.dmg` med Apple Developer ID-certifikat (`codesign`).

5. **Notarisering + stapling**
   - Skicka till Apple notarization (`notarytool`) och stapla ticket (`stapler`) så den fungerar utan Gatekeeper-varningar hos slutanvändare.

6. **Release/CI-pipeline**
   - Lägg till reproducerbar pipeline (GitHub Actions eller liknande) som bygger, testar, signerar, notariserar och bifogar `.dmg` till release.

## Kort rekommenderad väg

- Behåll `pyproject.toml` för Python-paket.
- Lägg till ett separat `scripts/build_macos_universal.sh` som:
  1) bygger `.app`,
  2) verifierar arkitekturer (`lipo -info`),
  3) signerar,
  4) notariserar,
  5) genererar `.dmg`.
- Kör detta i CI på macOS-runner med certifikat/hemligheter.
