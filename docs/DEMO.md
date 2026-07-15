# ReactorSim — demo-lathund

En sida att ha bredvid dig. Målet: 5–7 minuter som ser kontrollerat och proffsigt ut.

## Innan hon kommer (2 min)

1. **Stäng av Screen Time-nedtid.** ⚠️ Datorn har just nu en Screen Time-gräns
   (Downtime eller "Alla appar & kategorier") aktiv — den poppar en **"Time Limit"**-
   ruta på *varje* app med bundle, inklusive ReactorSim. Gå till
   **Systeminställningar → Skärmtid** och stäng av nedtiden/gränsen under demon.
   (Det är inget fel på appen — den nakna binären slapp bara för att Screen Time
   inte spårar oinpaketerade program.)
2. **Starta appen** från Spotlight: ⌘-mellanslag → "ReactorSim" → Enter. Ligger i
   `~/Applications/ReactorSim.app`. Öppnar på **PWR · Guided** vid full effekt.
3. **Ljud på**, volym uppåt halvvägs. Scram-tutan och PA-rösten är hela poängen.
4. Fönstret **öppnas maximerat** av sig självt nu. Vill du ändå ändra storlek går
   det med gröna knappen eller dubbelklick på titelraden (båda funkar numera).

## Wow-vägen (kör i den här ordningen)

| # | Gör | Säg |
|---|-----|-----|
| 1 | Låt henne titta. Peka på **reaktorkärlet** (RPV) till vänster — den lysande kapseln. | "Det där är själva härden. Färgen är den verkliga neutronflödesprofilen — varmast på mitten, kallare mot topp och botten." |
| 2 | Peka på siffrorna högst upp: 99,9 % effekt, 2 996 MWt, 988 MWe. | "Allt räknas ut från riktiga reaktorekvationer i realtid, inte en förinspelad animation." |
| 3 | Tryck **W** ett par gånger (dra ut stavarna lite). Titta på REACTOR POWER och REAKTIVITET. | "Drar jag ut styrstavarna stiger reaktiviteten och effekten — och kylningen svarar." |
| 4 | **SCRAM: tryck mellanslag.** Tutan tjuter, PA-rösten ropar "Reactor trip", stavarna faller, effekten kollapsar. | "Det här är en snabbstopp. Alla stavar in på två sekunder." |
| 5 | Vänta ~3 sek. Tryck **C** (kvittera larm), sedan **L** (återställ). | "Operatören kvitterar larmen och återställer — precis som i en riktig kontrollrum." |
| 6 | Tryck **U** (auto-startsekvens). Titta hur den kör upp anläggningen själv. Vill du visa lastföljning: slå på **O** (stav-auto) och dra ned **TBN GOV** (Q) — reaktorn följer ned och stabiliserar. | "Och den kan köra upp sig själv, och hålla T-avg automatiskt medan den följer turbinlasten." |
| 7 | Tryck **I** (haverimenyn) → välj t.ex. **DROPPED ROD** eller **STATION BLACKOUT**. | "Det finns en instruktörspanel med riktiga haverier — bortfallen stav, stationsblackout, ånggeneratorläcka…" |
| 8 | Tryck **M** ett par gånger — bläddra skinnen (Guided → Authentic Light → Authentic Dark). | "Samma anläggning, tre gränssnitt: det mjuka glasiga, och de platta som liknar riktiga ISA-101-HMI:er." |
| 9 | (Valfritt finale) Tryck **+** några gånger → ×600 tidskompression, se xenon och trender röra sig. | "Och jag kan snabba på tiden 600 gånger för att se de långsamma förloppen." |

## Om något krånglar

- **SCRAM-återställning gör inget?** Du tryckte nog mellanslag för snabbt. Ordningen är
  **mellanslag (scram) → C (kvittera) → L (reset)**. Reset är låst tills stavarna är
  helt inne *och* larmen kvitterade — det är meningen (så funkar en riktig anläggning).
- **Auto-start fastnar på "PRECLUDED — XENON"?** Xenonet efter ett scram kan förhindra
  kriticitet ett tag. Kör hellre auto-start (**U**) *innan* du scrammar, eller vänta ut det.
- **Tyst vid scram?** Kolla att ljudet är på i **Inställningar (,)** och att inte macOS är muteat.
- **Fel reaktor/skinn vid start?** Öppna **Inställningar (,)**, välj PWR + Guided.

## Fullständig tangentkarta

Stavar **W/S** · kylflöde **A/D** · turbinventil **Q/E** · matarvatten **F/V** ·
bor **B/G** · **mellanslag** SCRAM · **C** kvittera · **L** reset · **U** auto-start ·
**O** stav-auto · **I** haverier · **P** startillstånd · **T** turbintripp · **M** skinn ·
**,** inställningar · **+/−** hastighet · **Esc** paus · **F1–F6** dashboard-flikar.
Håll **Shift** för grova steg. Klicka på **kärlet** → härdkarta, på **NEUTRONICS**-rutan →
startpanel (1/M, ECP, stavvärde), på **brytarna** i ställverket för att manövrera dem.
