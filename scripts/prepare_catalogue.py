"""Extract the supplied catalogue without treating workbook text as instructions."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('workbook', nargs='?', type=Path,
                    default=ROOT / 'inputs' / 'hypercaloric_ganglia_datasets.xlsx')
SOURCE = parser.parse_args().workbook
(ROOT / 'metadata').mkdir(parents=True, exist_ok=True)
w = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
provenance = {'source_filename': SOURCE.name, 'sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'interpretation': 'User-supplied catalogue; external records and all sample claims require independent verification.'}
(ROOT / 'metadata' / 'catalogue_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
for sheet in w:
    with (ROOT / 'metadata' / f'catalogue_{sheet.title.lower().replace(" ", "_")}.csv').open('w', newline='') as f:
        csv.writer(f).writerows(sheet.iter_rows(values_only=True))
header, *rows = list(w['Datasets'].iter_rows(values_only=True))
records = [dict(zip(header, row)) for row in rows if row[0]]
with (ROOT / 'metadata' / 'dataset_audit.csv').open('w', newline='') as f:
    fields = ['id','accession','tissue','assay','candidate_role','verification_status','source_url']
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for r in records:
        acc = str(r['Accession'])
        writer.writerow(dict(id=r['ID'], accession=acc, tissue=r['Tissue (as deposited)'],
            assay=r['Assay'], candidate_role='diet comparison' if r['Tier']=='A' else 'context only' if r['Tier']=='B' else 'availability lead',
            verification_status='Not independently verified; NCBI access blocked at initial attempt',
            source_url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}' if acc.startswith('GSE') else ''))
print(f'Extracted {len(records)} catalogue records across {len(w.sheetnames)} sheets.')
