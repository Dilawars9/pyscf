# SPDX-License-Identifier: Apache-2.0
'''Assemble a self-contained offline HTML report from saved path analysis.'''

import argparse
import hashlib
import html
import json
from pathlib import Path

REPORT = Path(__file__).resolve().parents[1] / 'report'


def build_report(data, output):
    library = REPORT / 'vendor' / '3Dmol-2.5.5.min.js'
    license_text = (REPORT / 'vendor' / '3Dmol-LICENSE.txt').read_text()
    license_text += '\n\nBundled dependency notices:\n' + (REPORT / 'vendor' / '3Dmol-bundle-LICENSE.txt').read_text()
    license_text = '\n'.join(line.rstrip() for line in license_text.splitlines())
    # Prevent embedded text from terminating its HTML script element.
    encoded = json.dumps(data, allow_nan=False, separators=(',', ':')).replace('<', '\\u003c')
    viewer = library.read_text().replace('</script', '<\\/script')
    template = (REPORT / 'template.html').read_text()
    replacements = {'__REPORT_CSS__': (REPORT / 'report.css').read_text(),
                    '__VIEWER_LIBRARY__': viewer, '__DATA_JSON__': encoded,
                    '__REPORT_JS__': (REPORT / 'report.js').read_text(),
                    '__THIRD_PARTY_LICENSE__': html.escape(license_text),
                    '__VIEWER_SHA256__': hashlib.sha256(library.read_bytes()).hexdigest(),
                    '__IMAGE_COUNT__': str(data['metadata']['nimages'])}
    for marker, content in replacements.items():
        template = template.replace(marker, content)
    Path(output).write_text(template)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', type=Path)
    args = parser.parse_args()
    build_report(json.loads(args.analysis.read_text()), args.analysis.with_name('report.html'))
