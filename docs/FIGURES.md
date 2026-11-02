# Figure production specification

## Visual requirements agreed with the user

- Font: **Nimbus Sans** throughout. Italicize mouse gene symbols where appropriate.
- Control diet: **#0072B2**, RGB 0,114,178, blue circles.
- Hypercaloric diet: **#CC79A7**, RGB 204,121,167, pink triangles.
- Colors were read from the Okabe–Ito row on page 2 of the user-supplied
  `Colorblind Safe Color Schemes.pdf` (NCEAS, updated June 29, 2022).
  The reference PDF itself is not redistributed here.
- Distinguish groups with text/shape as well as color.
- Use 89-mm single-column or 183-mm double-column layouts as appropriate,
  7-pt text (6.5-pt ticks/legends), lowercase bold panel labels and 0.6-pt axes.
- Export vector PDF with valid embedded fonts, editable SVG and 600-dpi PNG.
  Check the specific journal's current requirements at submission.

## Planned data figures

1. **Primary diet-effect figure:** separate Pnoc and Oprl1 forest plots with
   log2(hypercaloric/control) effects, nominal 95% confidence intervals,
   biological n, and explicit ganglion, exposure/time and assay labels.
   No pooled effect across incompatible exposures or assays.
2. **Biological-replicate expression panels:** individual normalized expression
   values for each group and target. Identify the normalization and any plotting
   transform. Respect pooled biological units. Do not infer fold changes from
   arbitrary plotting pseudocounts.
3. **Single-cell supporting figure, if eligible:** replicate-aware cell-type
   localization/detection and pseudobulk diet effects; disclose cell counts,
   independent biological units, sparse targets and batch limitations.
4. **Supporting QC and secondary-gene panels:** sample similarity, library-size
   diagnostics and Pcsk1/Pcsk2/Cpe context where informative.

Each figure must have a source-data CSV/TSV and a full legend specifying n, unit
of replication, test, adjustment family, confidence interval and exclusions.
No biological data figure has been generated at handoff.
