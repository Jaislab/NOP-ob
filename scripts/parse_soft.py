"""Parse GEO family SOFT files into per-sample metadata tables (no interpretation)."""
from pathlib import Path
import gzip, sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def parse(acc):
    text = gzip.decompress((ROOT/'metadata'/acc/f'{acc}_family.soft.gz').read_bytes()).decode(errors='replace')
    rows, cur = [], None
    for line in text.splitlines():
        if line.startswith('^SAMPLE'):
            cur = {'gsm': line.split('=',1)[1].strip()}; rows.append(cur); continue
        if line.startswith('^'):
            cur = None; continue
        if cur is None or not line.startswith('!Sample_'): continue
        k, _, v = line[8:].partition(' = ')
        if k.startswith('characteristics_ch1') and ':' in v:
            kk, _, vv = v.partition(':'); k = 'char_' + kk.strip().lower().replace(' ','_')
        if k in cur: cur[k] = cur[k] + ' | ' + v
        else: cur[k] = v
    df = pd.DataFrame(rows)
    keep = [c for c in df.columns if c in ('gsm','title','source_name_ch1','description','supplementary_file_1','relation','library_strategy','molecule_ch1','extract_protocol_ch1','data_processing') or c.startswith('char_')]
    df[keep].to_csv(ROOT/'metadata'/acc/'samples_geo.csv', index=False)
    return df[keep]

if __name__ == '__main__':
    for acc in sys.argv[1:]:
        df = parse(acc)
        print('=====', acc, df.shape)
        with pd.option_context('display.max_columns', None, 'display.width', 250, 'display.max_colwidth', 60):
            print(df.drop(columns=[c for c in ['relation','extract_protocol_ch1','data_processing','supplementary_file_1'] if c in df]).to_string())
