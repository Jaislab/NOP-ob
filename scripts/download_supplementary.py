"""Download deposited GEO supplementary expression files with SHA-256 provenance."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os, sys, time
import requests

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    'GSE225224': ['GSE225224_RAW.tar', 'filelist.txt'],
    'GSE198895': ['GSE198895_RAW.tar', 'filelist.txt'],
    'GSE281091': ['GSE281091_Partek_Bulk_RNA-Seq_Exp_3A_Normalization_Normalized_counts_metabolic_states.txt.gz',
                  'GSE281091_Partek_Bulk_RNA-Seq_Exp_3B_Normalization_Normalized_counts_diet_types.txt.gz'],
    'GSE320302': ['GSE320302_rsem_count.txt.gz'],
    'GSE233163': ['GSE233163_cells_metadata.csv.gz', 'GSE233163_cells_counts.csv.gz'],
    'GSE309905': ['filelist.txt'],
}

def main(accs):
    s = requests.Session()
    s.verify = os.environ.get('REQUESTS_CA_BUNDLE') or os.environ.get('SSL_CERT_FILE') or True
    s.headers['User-Agent'] = 'NOP-diet-study/1.0 (public academic expression reanalysis)'
    log_path = ROOT/'metadata'/'download_manifest.json'
    log = json.loads(log_path.read_text()) if log_path.exists() else {}
    for acc in accs:
        for name in FILES[acc]:
            url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{name}'
            dest = ROOT/'data'/'raw'/acc/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            for attempt in range(5):
                try:
                    with s.get(url, stream=True, timeout=(20, 300)) as r:
                        r.raise_for_status()
                        h = hashlib.sha256(); n = 0
                        with open(dest, 'wb') as f:
                            for chunk in r.iter_content(1 << 20):
                                f.write(chunk); h.update(chunk); n += len(chunk)
                    break
                except requests.RequestException as e:
                    print('retry', name, e); time.sleep(2 ** (attempt + 1))
            else:
                raise SystemExit(f'failed {url}')
            log[f'{acc}/{name}'] = {'url': url, 'bytes': n, 'sha256': h.hexdigest(),
                                    'retrieved_utc': datetime.now(timezone.utc).isoformat()}
            print(acc, name, n)
            log_path.write_text(json.dumps(log, indent=2) + '\n')
            time.sleep(2)

if __name__ == '__main__':
    main(sys.argv[1:] or list(FILES))
