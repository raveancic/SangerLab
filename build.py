#!/usr/bin/env python3
"""Build the dependency-free, offline-capable SangerLab HTML from reviewed sources."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    source = ROOT / 'src'
    reference = ''.join(line.strip() for line in (ROOT / 'assets' / 'rCRS_15990.fasta').read_text().splitlines() if not line.startswith('>'))
    if len(reference) != 16569:
        raise ValueError('Expected the verified complete rCRS rotation.')
    content = (source / 'template.html').read_text()
    pieces = {
        'STYLE': (source / 'style.css').read_text(),
        'REFERENCE': json.dumps({'name': 'rCRS_15990', 'origin': 15990, 'sequence': reference}),
        'CORE': (source / 'core.js').read_text(),
        'ABIF': (source / 'abif.js').read_text(),
        'APP': (source / 'app.js').read_text(),
    }
    for name, value in pieces.items():
        content = content.replace('/*__' + name + '__*/', value)
    if '/*__' in content:
        raise ValueError('Unreplaced template placeholder.')
    license_text = (ROOT / 'LICENSE').read_text()
    content = content.replace('<html lang="it">', '<!--\n' + license_text + '\n-->\n<html lang="it">', 1)
    (ROOT / 'index.html').write_text(content)
    print(f'Built {ROOT / "index.html"} ({len(content.encode()):,} bytes)')


if __name__ == '__main__':
    main()
