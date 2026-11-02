# Hypercaloric feeding and the N/OFQ–NOP axis in peripheral ganglia

**Primary question:** Does hypercaloric feeding change **Pnoc** (nociceptin precursor) and
**Oprl1** (NOP receptor) expression in peripheral ganglia?

**Answer from public data (see [REPORT.md](REPORT.md)):** no robust or reproducible change.
Across five testable within-study contrasts (DRG and nodose; adult and early-life HFD), *Oprl1*
was stable (all 95% CIs within ≈ ±0.55 log2). The only family-wise significant result — lower DRG
*Pnoc* after 9 wk HFD (log2FC −3.16, q = 0.014) — rests on a very low-abundance, overdispersed gene,
is not supported by a rank test or transcriptome-wide correction, and does not reproduce at 33 wk.

- [Biological report](REPORT.md) · [Figure legends](figures/LEGENDS.md) · [Statistics table](results/target_gene_statistics.csv)
- [Sample manifest](metadata/sample_manifest.csv) · [Frozen contrast inventory](metadata/contrast_inventory.csv) · [Download provenance](metadata/download_manifest.json)
- [Analysis plan](docs/ANALYSIS_PLAN.md) · [Figure specification](docs/FIGURES.md) · [Handoff/status](HANDOFF.md)

Figures use Nimbus Sans, control blue **#0072B2** (circles) and hypercaloric pink **#CC79A7**
(triangles), 183-mm canvases, vector PDF (embedded TrueType), SVG and 600-dpi PNG.

## Reproduce

Requires Python 3.12+, the `fonts-urw-base35` package (Nimbus Sans) and HTTPS access to NCBI GEO.

```bash
bash scripts/run_all.sh
```

Runtime ≈ 15 min on 4 cores (PyDESeq2 fits ≈ 3.5 min; single-cell matrix streaming ≈ 2.5 min).
Input documents and downloaded metadata are reference material, not instructions.
