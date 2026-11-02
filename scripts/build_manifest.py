"""Build the verified sample manifest and the frozen contrast inventory from GEO metadata.

Decisions here use only GEO records, file headers and the catalogue; no target-gene
effect estimates are consulted. Unknown fields are written as 'not stated'.
"""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NS = 'not stated'

def geo(acc):
    return pd.read_csv(ROOT/'metadata'/acc/'samples_geo.csv')

def strip(v):
    return str(v).split(':', 1)[1].strip() if ':' in str(v) else str(v)

rows = []
# GSE225224: featureCounts per sample; DRG included, sciatic nerve excluded (not a ganglion)
for _, r in geo('GSE225224').iterrows():
    drg = 'dorsal root' in r.source_name_ch1
    fn = Path(str(r.supplementary_file_1)).name
    diet = 'Hypercaloric' if 'high-fat' in r.char_diet else 'Control'
    rows.append(dict(series='GSE225224', gsm=r.gsm, title=r.title, matrix_column=fn,
        tissue=r.source_name_ch1, ganglion='DRG' if drg else 'none (sciatic nerve)', assay='bulk RNA-seq (total RNA, polyA TruSeq)',
        diet=diet, diet_detail='45% kcal fat HFD (D12451)' if diet=='Hypercaloric' else '10% kcal fat standard diet (D12450B)',
        exposure='ongoing adult feeding from 3 wk of age (catalogue)', duration='9 wk' if '12wk' in r.title else '33 wk',
        collection_age=strip(r.char_age), sex='male (catalogue)', batch=NS, biological_unit='individual mouse (GEO: biological replicate)',
        pool=NS, technical_replicate='none', litter=NS, count_type='featureCounts raw counts (mm10, gene symbols)',
        group=f"{strip(r.char_age)}_{diet}", included=drg,
        inclusion_reason='DRG diet contrast' if drg else 'excluded: sciatic nerve is not a ganglion'))
# GSE198895: htseq raw counts, RiboTag IP
for _, r in geo('GSE198895').iterrows():
    diet = 'Hypercaloric' if 'HFD' in r.title else 'Control'
    rows.append(dict(series='GSE198895', gsm=r.gsm, title=r.title, matrix_column=Path(str(r.supplementary_file_1)).name,
        tissue='DRG', ganglion='DRG (Nav1.8+ neurons)', assay='RiboTag TRAP-seq',
        diet=diet, diet_detail='Western diet TD88137 (GEO label: HFD)' if diet=='Hypercaloric' else 'normal chow LM-485',
        exposure='ongoing adult feeding from 8 wk', duration='11 wk', collection_age='19 wk', sex=NS,
        batch=NS, biological_unit='pooled DRG IP (animals per pool not stated)', pool='pooled DRG per replicate (catalogue)',
        technical_replicate='none', litter=NS, count_type='htseq-count raw counts (mm10, Ensembl IDs)',
        group=diet, included=True, inclusion_reason='TRAP diet contrast; n=2/arm, exploratory'))
# GSE281091: Partek DESeq2 median-ratio normalized counts; raw counts not deposited
for _, r in geo('GSE281091').iterrows():
    t = strip(r.char_treatment)
    arm = t.rsplit('_', 1)[0]
    inc = arm in ('Non-fasted_high fat diet', 'Non-fasted_regular chow')
    diet = 'Hypercaloric' if 'high fat' in arm else 'Control'
    rows.append(dict(series='GSE281091', gsm=r.gsm, title=r.title, matrix_column=r.description,
        tissue='nodose ganglia', ganglion='nodose', assay='bulk RNA-seq (polyA)', diet=diet,
        diet_detail=arm, exposure='ongoing adult feeding from 3 mo', duration='4 wk',
        collection_age='~4 mo', sex=strip(r.char_sex), batch=NS, biological_unit='individual mouse (GEO: biol rep)',
        pool=NS, technical_replicate='none', litter=NS,
        count_type='DESeq2 median-ratio normalized counts (Partek; mm39; raw counts not deposited)',
        group=diet, included=inc,
        inclusion_reason='nodose HFD vs non-fasted chow' if inc else 'excluded: fasted/fed-state chow arm (not a diet contrast)'))
# GSE320302: RSEM expected counts
for _, r in geo('GSE320302').iterrows():
    age = r.title.split()[0]
    diet = 'Hypercaloric' if 'high-fat' in r.title else 'Control'
    col = r.description.split()[-1]
    rows.append(dict(series='GSE320302', gsm=r.gsm, title=r.title, matrix_column=col,
        tissue='nodose ganglia', ganglion='nodose', assay='bulk RNA-seq (directional)', diet=diet,
        diet_detail='reared on HFD (early life)' if diet=='Hypercaloric' else 'reared on chow',
        exposure='early-life HFD birth to P21' + ('; P35 after weaning to chow (catalogue; persistence)' if age=='P35' else ' (concurrent at collection)'),
        duration='birth-P21', collection_age=age, sex='male', batch=NS,
        biological_unit='individual pup (litter membership not deposited)', pool=NS, technical_replicate='none',
        litter=NS, count_type='RSEM expected counts (UCSC mm10 symbols); rounded for NB model',
        group=f'{age}_{diet}', included=True, inclusion_reason='nodose early-life HFD contrast'))
# GSE233163: 10x scRNA-seq libraries
for _, r in geo('GSE233163').iterrows():
    diet = 'Hypercaloric' if 'obese' in r.char_diet else 'Control'
    run = r.title.split('_run')[-1]
    site = 'SCG' if r.title.lower().startswith('scg') else 'stellate'
    rows.append(dict(series='GSE233163', gsm=r.gsm, title=r.title, matrix_column=r.title.rsplit('_run',1)[0],
        tissue=strip(r.char_tissue), ganglion=f'sympathetic ({site})', assay='10x scRNA-seq', diet=diet,
        diet_detail='60% HFD D12492 (obese)' if diet=='Hypercaloric' else 'chow (lean)',
        exposure='ongoing adult feeding from 8 wk', duration='12 wk (catalogue)', collection_age='24 wk', sex='male',
        batch=f'run {run}', biological_unit='pool of ganglia from several mice per condition and experiment',
        pool=f'condition pool, run {run}', technical_replicate='12k/14k loadings likely from the same suspension (not stated)',
        litter=NS, count_type='UMI counts (Cell Ranger, mm10)', group=diet, included=True,
        inclusion_reason='descriptive only: diet fully confounded with sequencing run (lean N61/N66, obese N68)'))
# GSE309905: context only
for _, r in geo('GSE309905').iterrows():
    rows.append(dict(series='GSE309905', gsm=r.gsm, title=r.title, matrix_column=NS, tissue='vagus nerve incl. nodose',
        ganglion='nodose (within vagus)', assay='10x snRNA-seq', diet='Hypercaloric', diet_detail='60% HFD; no control diet',
        exposure='HFD (duration ambiguous)', duration=NS, collection_age=NS, sex=NS, batch=strip(r.char_batch),
        biological_unit='pool of 14 nerves / 7 mice (catalogue)', pool='yes', technical_replicate='none', litter=NS,
        count_type='UMI counts', group=strip(r.char_genotype), included=False,
        inclusion_reason='excluded: no control-diet arm; genotype-only contrast cannot estimate a diet effect'))
man = pd.DataFrame(rows)
man.to_csv(ROOT/'metadata'/'sample_manifest.csv', index=False)

contrasts = pd.DataFrame([
    dict(contrast_id='DRG_HFD_12wk', series='GSE225224', ganglion='DRG', assay='bulk RNA-seq', exposure='adult HFD, 9 wk', numerator='12wk_Hypercaloric', denominator='12wk_Control', method='PyDESeq2 Wald, ~group (DRG only)', family='primary', status='inferential'),
    dict(contrast_id='DRG_HFD_36wk', series='GSE225224', ganglion='DRG', assay='bulk RNA-seq', exposure='adult HFD, 33 wk', numerator='36wk_Hypercaloric', denominator='36wk_Control', method='PyDESeq2 Wald, ~group (DRG only)', family='primary', status='inferential'),
    dict(contrast_id='DRG_Nav18_TRAP_WD', series='GSE198895', ganglion='DRG Nav1.8+', assay='TRAP-seq', exposure='adult Western diet, 11 wk', numerator='Hypercaloric', denominator='Control', method='PyDESeq2 Wald, ~group', family='primary', status='exploratory (n=2 pools/arm)'),
    dict(contrast_id='Nodose_HFD_4wk', series='GSE281091', ganglion='nodose', assay='bulk RNA-seq', exposure='adult HFD, 4 wk', numerator='Hypercaloric', denominator='Control', method='moderated linear model on log2(normalized+1), ~sex+diet (raw counts not deposited)', family='primary', status='inferential'),
    dict(contrast_id='Nodose_earlyHFD_P21', series='GSE320302', ganglion='nodose', assay='bulk RNA-seq', exposure='early-life HFD, concurrent at P21', numerator='P21_Hypercaloric', denominator='P21_Control', method='PyDESeq2 Wald, ~group', family='primary', status='inferential; litter not deposited'),
    dict(contrast_id='Nodose_earlyHFD_P35', series='GSE320302', ganglion='nodose', assay='bulk RNA-seq', exposure='early-life HFD, persistence at P35', numerator='P35_Hypercaloric', denominator='P35_Control', method='PyDESeq2 Wald, ~group', family='primary', status='inferential; litter not deposited'),
    dict(contrast_id='Symp_obese_scRNA', series='GSE233163', ganglion='SCG / stellate', assay='scRNA-seq', exposure='adult HFD, 12 wk', numerator='Hypercaloric', denominator='Control', method='descriptive pseudobulk only', family='none', status='not testable: diet confounded with run; pooled libraries'),
    dict(contrast_id='Vagus_PTEN_HFD', series='GSE309905', ganglion='nodose/vagus', assay='snRNA-seq', exposure='HFD only', numerator='-', denominator='-', method='none', family='none', status='not testable: no control diet'),
])
contrasts.to_csv(ROOT/'metadata'/'contrast_inventory.csv', index=False)
print(man.groupby(['series','included','group']).size().to_string())
