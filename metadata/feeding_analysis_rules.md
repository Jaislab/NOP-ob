# Feeding-regulation analysis (modules A-C): prespecified rules

Frozen before any fasting/refeeding or panel effect was estimated (the Pnoc/Oprl1 HFD results
from the primary analysis were already known).

- GSE281091: the two deposited normalized files (3A: fasted, refed, non-fasted; 3B: HFD, non-fasted)
  share the six non-fasted samples. 3A values are rescaled onto the 3B scale with a per-sample
  factor from the shared samples (median ratio over genes present in both); analysis uses genes
  present in both files. One moderated (limma-trend-style) linear model, ~group + sex, 24 samples.
- BH families: A_nop (Pnoc, Oprl1 x fasted-vs-fed, refed-vs-fasted); A_panel (vagal_feeding x the
  same two contrasts); A_ieg (ieg x refed-vs-fasted); C_panel (vagal_feeding x HFD adult, P21, P35).
  Panel genes absent or with mean normalized count < 5 are reported as not testable and excluded
  from the family.
- Literature-predicted directions (rat nodose) are recorded in feeding_gene_panels.csv and are
  evaluated descriptively (sign agreement) alongside the two-sided tests.
- Module B: see feeding_contrast_inventory.csv. Known risk: the three HFD-male samples sit on a
  technical PC1 axis; the female-only sensitivity is reported regardless of outcome.
