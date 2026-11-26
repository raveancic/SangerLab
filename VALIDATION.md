# Verifica della versione pilota 0.1.0

Data: 26 novembre 2026. Questi risultati descrivono i controlli eseguiti su questa
versione; non costituiscono una validazione generale per qualsiasi strumento ABI,
qualsiasi campione o qualsiasi tablet.

## Motore e parser

- 31 test JavaScript superati: 16 per il parser ABIF e 15 per il motore.
- Il parser è stato confrontato con Biopython su un cromatogramma reale: basi,
  Phred, posizioni dei picchi e ogni valore delle quattro tracce coincidono.
- 49 casi di allineamento confrontati con `Bio.Align.PairwiseAligner`:
  punteggi affini, traceback, orientamento, coordinate, qualità e alleli.
- 132 casi di trimming confrontati con un algoritmo indipendente che enumera
  tutti gli intervalli possibili delle sequenze sintetiche.
- Verificati entrambi i punti di giunzione: origine del FASTA ruotato e origine
  della numerazione canonica rCRS. Le due giunzioni sono concetti distinti.
- Verificati l'invarianza delle coordinate SNP dopo reverse complement e cambio
  della rotazione, gli indel nelle ripetizioni e l'assenza di Phred inventati per
  delezioni o tag PCON assenti.

Ambiente dei controlli: Node.js 24.22.3, Python 3.12, Biopython 1.88.
I test riproducibili e le istruzioni sono in `tests/` e `README.md`.
La versione standalone del programma non richiede queste dipendenze.

## Interfaccia: cosa è stato e cosa non è stato provato

Sono stati eseguiti 15 controlli di flusso in un DOM simulato (JSDOM), con
dimensioni di canvas di 1.100 e 330 pixel. I controlli comprendono caricamento
del file, trimming, chiamata del motore mediante un adattatore Worker per Node,
rimando da tabella a picco, numerazione all'origine, tab dell'allineamento,
cinque esportazioni, invalidazione dei risultati dopo modifiche e rifiuto di
reference incompatibili. Tutti sono stati superati. Il disegno della traccia è
stato anche ispezionato come raster prodotto dal medesimo codice Canvas.

**Questa simulazione non è un test in un browser reale.** Non verifica il layout
CSS, l'applicazione della Content Security Policy, l'avvio di un Worker nativo,
i selettori file o i download dei browser mobili, né il touch. L'accesso a un
browser per il collaudo grafico non era disponibile in questo ambiente.

## Collaudo da completare prima dell'aula

1. Aprire `index.html` in Chrome o Firefox desktop e caricare un file reale.
2. Pubblicare il programma su HTTPS e aprirlo sul tablet Android del corso.
3. Verificare importazione dal gestore file, lettura dei Phred, selezione delle
   basi con il dito, frecce di navigazione e slider.
4. Confrontare coordinate e differenze con il risultato indipendente atteso.
5. Scaricare FASTA, TSV e JSON e riaprirli dal gestore file del tablet.
6. Ripetere su almeno un secondo dispositivo e, se previsto, Safari su iPad.

Usare letture forward/reverse e almeno un campione con estremità scadenti,
un'indel in una ripetizione e una sequenza non omologa. Questa versione elabora
singole letture: l'assemblaggio bidirezionale è una funzionalità futura.

I dati individuali e i relativi report non fanno parte di questo repository.

## Aggiornamento cumulativo 0.1.1

La versione pubblica 0.1.1 integra tutte le revisioni locali dell'interfaccia
(Neo-Brutalism, confronto reference/campione e aiuto Mott/Phred). Il motore di
analisi e il parser non cambiano, salvo il numero di versione del report.

- 31/31 test Node.js superati nuovamente (Node.js 22.16.0).
- Chrome/Chromium headless: documento caricato direttamente nel browser tramite
  Playwright; entrambi i dialoghi si aprono e si chiudono, il file reale
  Example7.ab1 conserva 857 basi con p=0,01, orientamento forward, coordinate
  16055–341 e 11 differenze candidate.
- Secondo test browser con un AB1 sintetico costruito come reverse complement
  del file reale, incluse qualità, coordinate dei picchi e canali invertiti e
  complementati. Stessi risultati canonici; verificate etichette reference e
  campione, coordinate AB1 decrescenti, selezione via click del picco C16069T,
  navigazione e invalidazione dell'allineamento cambiando il cutoff a 0,05
  (903 basi conservate). Nessuna eccezione JavaScript nel test reverse.
- Ispezionato il raster del canvas reverse prodotto dal browser.

Limiti: il browser dell'ambiente blocca la navigazione diretta a URL file://;
le prove sono state effettuate inserendo il documento HTML in Chromium, non
visitando il sito GitHub Pages. La compatibilità su tablet Android/iPad fisici
resta da verificare. Il cromatogramma reale e la sua derivazione sintetica non
sono distribuiti nell'aggiornamento.

## v0.1.2: visibilita della reference nel cromatogramma

31/31 test Node superati. Chromium headless: HTML caricato come contenuto
(in questo ambiente la navigazione diretta file:// e bloccata), AB1 reale
Example7.ab1, trimming 857 basi, coordinate 16055-341, 11 differenze.
Controllate le istruzioni di disegno Canvas dopo allineamento e il pulsante
Vedi picco: la C della reference e la T della lettura condividono la stessa
ascissa nella variante C16069T. Screenshot verificato visivamente.
Non testato su dispositivo Android fisico o sul sito GitHub Pages aggiornato.

## v0.1.3 — Contrasto, tema scuro e loghi

- 31 test Node del parser e del core superati; del core cambia esclusivamente VERSION.
- Test Chromium headless con HTML caricato in una pagina di prova: Example7.ab1, 857 basi conservate, 11 differenze candidate; oggetto analysis identico prima e dopo il cambio di tema. Apertura popup di trimming e selezione del logo verificate.
- Viewport mobile 390 x 844 e input touch simulato: tema e selezione logo funzionanti, nessun overflow orizzontale del documento prima dell'analisi.
- Screenshot desktop Light/Dark e cromatogramma scuro ispezionati. La navigazione file:// e verso un dominio di prova e' bloccata dalle politiche del browser dell'ambiente; persistenza dopo ricaricamento, sito GitHub Pages e dispositivi Android/iPad fisici non verificati. La mancata disponibilita' di localStorage non impedisce l'uso dell'applicazione.

### v0.1.3 — grafica definitiva

Integrata esattamente l'anteprima approvata `SangerLab_design_logo_giallo.html`:
logo vettoriale fisso, giallo #FFD400 nei due temi, sfondo geometrico fisso
con quattro forme, pannelli bianchi o antracite. Nessun selettore del logo.
Verificata l'identità dei quattro blocchi script/reference con l'anteprima
precedente e l'identità byte per byte dell'HTML ricostruito con `build.py`
rispetto all'anteprima approvata. Controllo sintattico Node e 31/31 test
unitari superati. Il collaudo finale sul sito pubblicato e sul tablet
fisico resta da effettuare; nessun nuovo collaudo grafico in questa fase
che riguarda soltanto la preparazione del pacchetto.

## v0.1.4 — Marcatori delle differenze

- `node --test tests/*.test.cjs`: 31/31 superati; core e parser invariati.
- Chromium, HTML caricato in memoria: Example7.ab1 produce 11 differenze candidate; verificato un asterisco dopo «Vedi picco» sulla prima differenza, mantenimento del simbolo al cambio tema, assenza nella regione non allineata e cancellazione dopo reset del trimming.
- Controllo visivo dello screenshot in Dark mode: asterisco sotto la base campione corrispondente, senza sovrapposizione alla traccia.
- Non ripetuti in questa revisione i test grafici con letture reverse e delezioni reali, né il collaudo su tablet fisico. L'apertura diretta file:// è bloccata nell'ambiente di verifica; il test Chromium ha usato l'HTML in memoria.

## v0.1.5 — Navigazione diretta sul cromatogramma

- Derivata dalla v0.1.4; `src/core.js` e `src/abif.js` identici byte per byte. Versione applicativa aggiornata nel report JSON, versione del motore invariata.
- `node --check src/app.js` e `node --test tests/*.test.cjs`: 31/31 superati.
- Chromium con HTML caricato in memoria e Example7.ab1: allineamento accettato, 11 differenze candidate; verificati eventi wheel orizzontali in entrambi i versi, trascinamento mouse, clic di ispezione, wheel verticale senza spostamento della sequenza, swipe orizzontale touch simulato e tap di ispezione. Il risultato dell'allineamento e' identico prima e dopo le interazioni.
- Trackpad fisico, tablet fisico e navigazione diretta con una lettura reverse reale non ancora collaudati. Il calcolo dello spostamento usa il verso della traccia visualizzata. Il comportamento nativo della barra inferiore resta invariato rispetto alla v0.1.4.

## v1.0.6

Presentation-only update: compact English README, icon-only theme switch, removed visible version and local badge. Scientific core and ABIF parser unchanged. All 31 Node tests passed.
