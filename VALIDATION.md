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
