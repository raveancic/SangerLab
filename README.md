# SangerLab v1.0.6

SangerLab is a lightweight, browser-based tool for Sanger sequencing analysis, developed for the Conservation Genetics course at the University of Pavia.

The interface is currently in Italian, as the software is designed for practical classes taught in Italian.

## Try it online

https://raveancic.github.io/SangerLab/

## Features

- AB1 file import and interactive chromatogram visualization.
- Phred quality scores and automatic/manual trimming.
- Automatic forward/reverse complement detection.
- Circular mitochondrial DNA alignment against the rCRS (NC_012920.1).
- Canonical rCRS coordinate mapping.
- Identification and visualization of candidate substitutions and indels.
- Interactive navigation using mouse, trackpad or touchscreen.
- FASTA, TSV and JSON export.
- Light and dark modes.

## Usage

Open the web application or run `index.html` locally in a modern browser. No installation, registration or external dependencies are required.

All analyses are performed locally in the browser. Sequence data are never uploaded to external servers.

The current version supports single-read analysis of human mitochondrial DNA. Paired forward/reverse read assembly and consensus generation are planned for future releases.

## License

MIT License.

The AB1 parser is adapted from [Sanger Viewer](https://github.com/animeshkundu/sanger-viewer).
