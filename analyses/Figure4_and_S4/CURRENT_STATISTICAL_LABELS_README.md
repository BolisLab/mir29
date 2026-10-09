# Current statistical annotations — 9 October 2026

`CURRENT_STATISTICAL_LABELS.tsv` maps the verified statistical outputs to current
Figure 4/S4 annotations and reference-distance support. It distinguishes:

- 4d pooled / S4a matched global treatment-by-categorical-time interactions.
- 4f pooled T6 entropy/FAE contrasts and matched alternatives reported in its
  legend; these are endpoint comparisons, not global curve tests.
- 4i/j pooled / S4c/e matched global average-shift tests.
- S4d matched global average shifts, BH across two arms separately per clock.
- S4k matched hTERT RNA-clock global ML likelihood-ratio shift tests, retaining
  BH across two shifts plus two interactions per clock.
- 4o centre: matched 123-CpG reference-distance arithmetic mean T1–T6;
  4o right: matched mean T5–T6. Both use three run summaries per arm and BH
  across two arms separately for each endpoint. q = 0.025 is not a slope test.
- hTERT reference support: the same 123 CpGs/medians at the third-cycle 96-h
  endpoint and, separately, mean over all three cycles. These are not additional
  panels in Figure S4 and must not be labelled primary-HGF T5–T6 tests.

Figure 4l is unchanged. Clock curves display global q values only, without
individual-timepoint stars. Historical post-hoc tables remain available but do
not determine current graphical annotations.

Sources are `tables/main_statistics_own_vs_pooled_long.tsv`,
`tables/S4_epigenetic_clock_global_statistics.tsv`,
`tables/S4k_clock_global_statistics.tsv`, `tables/4o_global_tests.tsv`,
`tables/4o_late_tests.tsv` and `tables/hTERT_123CpG_reference_tests.tsv`.
Executable replay is in code modules 03, 05b and 06a–06c; annotation verification
and display are in module 08. No new test or correction family is introduced.
