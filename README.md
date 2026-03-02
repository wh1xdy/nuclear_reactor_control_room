# Nuclear Reactor Control Room Sim

En kontrollrumssimulator i Python/Pygame med tre reaktormodeller:

- PWR
- BWR
- RBMK

Nu innehåller projektet också ett **supervisor-lager** som lägger till:

- förenklad balance-of-plant (tryck, ånginventarie, kondensor, matarvatten),
- skyddssystem med I&C-kanaler (2oo3-votering, filter/fel), alarmer och latched SCRAM,
- startinterlock (startup permit),
- injicerbara fel (pumpdegradering, matarvattenförlust),
- axial tvåfas-kanalmodell i BWR/RBMK (kvalitet, slip-baserad void och tryckfall över kanal).

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
- `F/V`: matarvattenventil
- `H/N`: pressurizer heater
- `P`: startup permit interlock (på/av)
- `T`: turbine trip (på/av)
- `Z`: pumpfel (på/av)
- `X`: matarvattenfel (på/av)
- `C`: kvittera aktiva larm
- `L`: återställ SCRAM-latch (permissive-matris)
- `B/G`: bypass för högtrycks-/högbränsletrip
- `I`: inhibit auto-SCRAM
- `SHIFT`: snabb justering
- `SPACE`: SCRAM
- `R`: reset
- `ESC`: avsluta

## Validering

Kör kvalitetssäkring av transienta förlopp:

```bash
python validation.py
```

Skriptet verifierar kvalitativt att:

1. SCRAM minskar effekt i alla tre reaktortyper.
2. Skyddssystem och trips fungerar deterministiskt.
3. BWR/RBMK-transienter körs med axial tvåfas-void/tryckfallsmodell.
4. Referenstransienter jämförs mot kalibreringskurvor (`data/reference_transients.json`).

## Realism och begränsningar

Detta är fortfarande en utbildningssimulator, inte en certifierad träningssimulator.
Men jämfört med MVP-versionen ingår nu fler realistiska kontroll- och skyddskedjor
samt fler sekundärsystemsinteraktioner.

## macOS universal `.dmg`

Repo:t har CI-flöde som bygger två `.app`-bundles (x86_64 + arm64), slår ihop dem till
en universal `.app` och paketerar en DMG-artefakt.

Kvar för en full distributionskedja för slutanvändare:

1. **Kodsignering** av `.app` och `.dmg` (`codesign`, Developer ID).
2. **Notarisering** + stapling (`notarytool` + `stapler`).
3. **Releasepublicering** (t.ex. bifoga DMG till GitHub Release).
