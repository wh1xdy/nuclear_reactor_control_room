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
