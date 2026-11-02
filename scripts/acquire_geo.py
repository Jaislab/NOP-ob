"""Retrieve official GEO metadata; preserve TLS verification and inherited proxy.

This only acquires records. Sample group assignments and file semantics must be
reviewed before any differential-expression model is fitted.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import argparse
import gzip
import hashlib
import json
import os
import re
import requests

ROOT = Path(__file__).resolve().parents[1]
ACCESSIONS = ['GSE225224','GSE198895','GSE281091','GSE320302','GSE233163','GSE309905']

def fetch(acc):
    if not re.fullmatch(r'GSE\d+', acc):
        raise ValueError(f'Invalid accession: {acc}')
    prefix = acc[:-3] + 'nnn'
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}'
    dest = ROOT / 'metadata' / acc
    dest.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    # Use the environment's CA bundle if supplied. Never disable verification.
    session.verify = os.environ.get('REQUESTS_CA_BUNDLE') or os.environ.get('SSL_CERT_FILE') or True
    session.headers['User-Agent'] = 'NOP-diet-study/1.0 (public academic expression reanalysis)'
    outcomes = []
    for suffix, name in [(f'soft/{acc}_family.soft.gz', f'{acc}_family.soft.gz'),
                         ('suppl/', 'supplementary_listing.html')]:
        url = base + '/' + suffix
        path = dest / name
        entry = {'accession':acc, 'url':url, 'retrieved_utc':datetime.now(timezone.utc).isoformat()}
        try:
            response = session.get(url, timeout=(15, 90))
            response.raise_for_status()
            payload = response.content
            if name.endswith('.gz'):
                decoded = gzip.decompress(payload)
                if f'^SERIES = {acc}'.encode() not in decoded:
                    raise ValueError('Downloaded SOFT does not contain the requested series record')
            elif b'<html' not in payload.lower() and b'<!doctype' not in payload.lower():
                raise ValueError('Supplementary directory response is not HTML')
            path.write_bytes(payload)
            entry.update(status='downloaded', bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(), path=str(path.relative_to(ROOT)))
        except (requests.RequestException, ValueError, OSError) as exc:
            entry.update(status='failed', error=str(exc))
        outcomes.append(entry)
    return outcomes

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('accessions', nargs='*', default=ACCESSIONS)
    args = parser.parse_args()
    with ThreadPoolExecutor(max_workers=3) as executor:
        results = [item for group in executor.map(fetch, args.accessions) for item in group]
    path = ROOT / 'logs' / 'geo_acquisition.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results, indent=2) + '\n')
    for r in results:
        print(r['accession'], r['status'], r['url'])
    raise SystemExit(1 if any(r['status']=='failed' for r in results) else 0)

if __name__ == '__main__':
    main()
