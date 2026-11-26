# SangerLab — D-loop e coordinate rCRS

Piccolo strumento didattico per analizzare **una lettura Sanger mitocondriale alla volta**: cromatogramma AB1, qualità Phred, trimming e allineamento a una reference circolare, con numerazione rCRS originale.

Il file `index.html` contiene l'applicazione completa. Il calcolo avviene nel browser, senza account, backend, CDN o installazione di librerie. I file selezionati restano nel dispositivo: l'applicazione non li carica su un server.

## Avvio rapido

1. Apri `index.html` in un browser moderno, oppure apri il sito GitHub Pages pubblicato dal docente.
2. Premi **Apri .ab1** e seleziona un cromatogramma dal dispositivo.
3. La reference inclusa è la rCRS completa **NC_012920.1**, ruotata per iniziare dalla posizione canonica **15.990**. Con **Cambia FASTA** puoi usare la stessa rCRS in un'altra rotazione, inclusa quella canonica.
4. Controlla il trimming proposto. Puoi cambiare il cutoff Mott, modificare gli estremi manualmente oppure conservare tutte le basi.
5. Premi **Allinea alla rCRS**. La modalità automatica confronta la lettura e il suo reverse complement.
6. Esamina le differenze e usa **Vedi picco** per controllare il cromatogramma nella posizione corrispondente.
7. Esporta i risultati necessari prima di chiudere o ricaricare la pagina.

La sessione rimane in memoria e non viene salvata automaticamente.

### Uso offline e tablet

Il file HTML scaricato è autosufficiente: non scarica script, font o reference durante l'analisi. Su computer è possibile aprirlo direttamente dal filesystem. Alcuni browser o gestori di file mobili limitano l'apertura degli HTML locali o l'esecuzione di Worker da URL `blob:`; in questo caso usa la versione HTTPS su GitHub Pages.

Per un uso offline prevedibile, conserva una copia locale di `index.html`; la semplice visita al sito non installa una PWA e non garantisce che una futura ricarica senza rete funzioni.

Servono un browser moderno con JavaScript, File API, Canvas, Pointer Events, Blob, Web Workers e ResizeObserver. L'interfaccia offre pulsanti, selezione diretta delle basi, scorrimento e numero di basi visibili regolabile. **L'usabilità su un tablet fisico Android o iPad deve ancora essere verificata prima del laboratorio**: una prova su viewport mobile o con emulazione touch non la sostituisce.

## Quali analisi esegue

| Passaggio | Comportamento |
|---|---|
| Lettura AB1 | Legge basi, qualità, coordinate dei picchi e quattro tracce analizzate; usa l'ordine dei canali dichiarato nel file. |
| Qualità | Visualizza i Phred del basecaller presenti nel file; non li ricalcola dalla fluorescenza. |
| Cromatogramma | Mantiene il verso originale del file anche quando l'allineamento sceglie il reverse complement. |
| Trimming | Propone un intervallo Mott e permette di cambiare manualmente la prima e l'ultima base conservata. |
| Reference | Accetta una singola sequenza FASTA, completa, che sia una rotazione esatta della rCRS inclusa. |
| Allineamento | Cerca entrambe le orientazioni con allineamento semiglobale e gap affini su reference circolare. |
| Coordinate | Mostra posizione nel file AB1, posizione nella reference FASTA e posizione rCRS canonica. |
| Differenze | Riporta sostituzioni, indel e basi ambigue, con qualità e rimando al cromatogramma. |
| Esportazione | Lettura FASTA, allineamento FASTA, differenze TSV, qualità TSV e report JSON. |

Questa versione non modifica le basi, non richiama i picchi, non assembla letture forward/reverse, non genera un consenso di più letture, non esegue allineamenti multipli, non stima l'eteroplasmia e non assegna aplogruppi.

## Metodo e convenzioni scientifiche

### AB1 e Phred

Il parser legge preferenzialmente `PBAS2`, `PLOC2` e `PCON2`, cioè le chiamate originali del basecaller e i relativi picchi e qualità. Se la coppia originale di basi e picchi è completamente assente, può utilizzare la generazione 1 segnalando che può contenere modifiche precedenti. Le qualità di generazioni diverse non vengono mescolate.

I canali analizzati sono `DATA9`–`DATA12`; `FWO_1` ne definisce l'ordine. Tag incoerenti, lunghezze incompatibili o coordinate dei picchi fuori intervallo generano un errore. In assenza di Phred il programma segnala il limite e non inventa qualità; le analisi che li richiedono rimangono indisponibili.

Per un Phred `Q`, la probabilità nominale di errore della chiamata è `10^(-Q/10)`. La qualità di una singola base non è una probabilità della variante biologica e non dimostra da sola l'assenza di contaminazione o artefatti.

### Trimming Mott

Per ogni base viene calcolato:

```text
score(i) = cutoff − 10^(−Q(i)/10)
```

Si conserva il segmento contiguo con somma positiva massima. I cutoff disponibili sono **0,01**, **0,05** e **0,001**; il valore iniziale è **0,01**. Se nessun segmento ha somma positiva, il trimming automatico non propone una regione conservata.

Questo criterio rimuove estremità sfavorevoli, ma può conservare singole basi di qualità bassa all'interno di un segmento complessivamente buono. Non equivale a richiedere `Q ≥ 20` per ogni base conservata. Le qualità originali restano consultabili anche dopo il trimming.

Gli estremi mostrati nell'interfaccia sono **1-based, inclusivi**. Nel report JSON `trim.start` è invece **0-based inclusivo** e `trim.end` **0-based esclusivo**, come specificato anche nel report.

### Reference circolare e numerazione

La reference incorporata contiene **16.569 posizioni**, inclusa la `N` convenzionale alla posizione rCRS 3107. È una semplice rotazione di NC_012920.1: la sequenza biologica non cambia.

Se `s` è la posizione canonica della prima base della reference ruotata e `p` è una posizione FASTA 1-based, allora:

```text
posizione_rCRS = ((s + p − 2) mod 16569) + 1
```

Per la rotazione inclusa, `s = 15990`: le posizioni FASTA 580 e 581 corrispondono rispettivamente alle posizioni rCRS 16569 e 1. L'origine di un nuovo FASTA viene riconosciuta confrontandone l'intera sequenza con la rCRS, senza affidarsi all'intestazione.

Per consentire allineamenti che attraversano il punto di chiusura del FASTA, il motore usa due copie consecutive della reference. Riporta poi le coordinate modulo la lunghezza originale, conservando anche la distinzione tra attraversamento dell'origine FASTA e attraversamento dell'origine canonica rCRS.

FASTA frammentari, sequenze modificate mediante sostituzioni/indel o reference di altre specie vengono rifiutati da questa versione. Una semplice rotazione non richiede una trasformazione generale delle coordinate; una reference realmente alterata sì, e quel caso non è implementato.

### Allineamento e interpretazione delle differenze

La lettura conservata viene allineata per intero a un tratto della reference. Le regioni della reference prima e dopo l'allineamento non sono penalizzate; i gap all'interno della lettura lo sono. Punteggi:

| Elemento | Punteggio |
|---|---:|
| Identità fra basi A/C/G/T | +2 |
| Sostituzione fra basi A/C/G/T | −3 |
| Apertura di un gap | −6 |
| Ogni base aggiuntiva nello stesso gap | −1 |
| Confronto con simboli IUPAC ambigui | 0 |

La modalità automatica confronta entrambe le orientazioni. I criteri iniziali di controllo richiedono almeno 50 basi e identità di almeno 80%; sono filtri operativi, non una stima probabilistica di affidabilità. La tabella contiene **differenze candidate da controllare sui picchi**. Un allineamento segnalato come non accettato non diventa valido perché è esportabile.

Le inserzioni sono localizzate fra due basi di riferimento; non ricevono una nuova posizione rCRS intera. Le delezioni non hanno un Phred proprio: quando disponibile, viene riportata separatamente la qualità dei fianchi osservati. In ripetizioni o omopolimeri possono esistere posizionamenti equivalenti degli indel. Il programma segnala alcune ambiguità locali, ma **non applica normalizzazione mitocondriale, HGVS o VCF** e non enumera tutti gli allineamenti ottimali.

## File esportati

| Pulsante | Contenuto |
|---|---|
| Lettura FASTA | La lettura rifilata nell'orientamento della reference. Le regioni non osservate non vengono riempite con basi della reference. |
| Allineamento FASTA | Due record con gap: tratto della reference e lettura allineata. |
| Differenze TSV | Alleli, coordinate, qualità, ancore degli indel e indicatori di ambiguità. |
| Qualità TSV | Tutte le basi originali con Phred, indicazione conservata/esclusa e coordinate disponibili dopo l'allineamento. |
| Report JSON | Parametri, estremi del trimming, reference, risultati e metadati necessari a ricostruire l'analisi; comprende anche nome e hash del file, quando disponibili. |

L'output FASTA è una lettura singola curata mediante trimming, non un consenso forward/reverse. Gli output e i nomi dei campioni appartengono ai dati del laboratorio: restano separati dal repository pubblico del software.

## Pubblicazione gratuita su GitHub Pages

Per la modalità GitHub Free usa un repository pubblico. Pubblica **i file del programma**, mantenendo i file AB1 e i risultati del laboratorio sui dispositivi o nello spazio dedicato ai dati del corso.

1. Crea un repository, ad esempio `sangerlab`.
2. Estrai il pacchetto ZIP del progetto. Con **Add file → Upload files** carica il suo contenuto, mantenendo cartelle e file; `index.html` deve trovarsi nella radice del repository, non dentro un'ulteriore cartella e non soltanto dentro lo ZIP.
3. Conferma il caricamento sul branch `main`.
4. Apri **Settings → Pages**. In **Build and deployment**, imposta **Source → Deploy from a branch**.
5. Scegli **main** e **/(root)**, quindi **Save**.
6. Attendi la pubblicazione e apri il collegamento mostrato da GitHub. L'indirizzo avrà normalmente la forma `https://<utente>.github.io/sangerlab/`.

Gli studenti aprono il collegamento e selezionano i propri file nel browser; non devono avere un account GitHub. La pagina pubblicata distribuisce l'applicazione statica. Non è necessario eseguire Python o Node sul server.

I passaggi seguono la [documentazione GitHub per creare un sito Pages](https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site) e quella sulla [pubblicazione da un branch](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

## Sorgenti, build e test

`index.html` è già pronto. Il sorgente è suddiviso in `src/abif.js` (lettura AB1), `src/core.js` (algoritmi), `src/app.js` (interfaccia), `src/template.html` e `src/style.css`. La reference pubblica è in `assets/rCRS_15990.fasta`.

Per rigenerare l'HTML dopo una modifica ai sorgenti, dalla radice del progetto:

```bash
python3 build.py
```

Lo script usa soltanto la libreria standard Python. Per eseguire i test JavaScript con un Node.js moderno dotato di test runner:

```bash
node --test tests/abif.test.cjs tests/core.test.cjs
```

I test nel repository usano dati sintetici. Coprono parsing, qualità assenti, reverse complement, coordinate circolari, allineamento, trimming e input non validi. Il confronto indipendente con Biopython richiede un ambiente Python con Biopython installato e Node.js:

```bash
python3 tests/validate_biopython.py
```

Il validatore supporta anche dati reali locali, senza includerli nel repository. Consulta `python3 tests/validate_biopython.py --help` per gli argomenti. Conserva eventuali report prodotti con dati reali insieme ai dati stessi.

Il superamento dei test algoritmici non sostituisce il collaudo didattico: prima dell'uso in aula verifica importazione, selezione touch delle basi, navigazione, download dei risultati e prestazioni sui tablet effettivamente disponibili. Questa prima versione non è stata validata su tutti i modelli di sequenziatore o su tutte le configurazioni mobili.

## Licenza e fonti

Il codice è distribuito con licenza MIT; vedi `LICENSE` e `THIRD_PARTY_NOTICES.txt`.

Il parser è adattato da [Sanger Viewer](https://github.com/animeshkundu/sanger-viewer), revisione `3ceda33e18dde1700f6a53414f58ad4b1d8989be`, con controlli aggiuntivi. Trimming probabilistico, motore di allineamento e gestione delle coordinate sono implementati separatamente; non sono una copia dei rispettivi moduli upstream.

Fonti tecniche:

- [Applied Biosystems, ABIF File Format Specification and Sample File Schema, luglio 2006](https://archive.gfjc.fiu.edu/workshops/resources/literature/ABIF_File_Format.pdf).
- [NCBI RefSeq NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1).
- [Biopython — Bio.SeqIO.AbiIO](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py), usato come controllo indipendente del parsing.
- [Biopython — allineamento a coppie](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html), usato come controllo indipendente del punteggio di allineamento.

Biopython è una dipendenza facoltativa dei test indipendenti; non viene scaricato né eseguito dall'applicazione web.
