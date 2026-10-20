# Number versus combination of endorsed value domains and well-being: a cross-sectional survey of Japanese households affected by cancer and disaster

This repository holds the analysis code for the article "Number versus combination of endorsed value domains and well-being: a cross-sectional survey of Japanese households affected by cancer and disaster". It contains the code and the codebook of the analysis data file. The data file itself, the fitted model objects and the questionnaire are not included.

## Data

The scripts read the analysis data file at

    data/data_release_v53.csv

It has one row for each of the 2179 respondents and 92 columns: a serial number, the covariates, the K6 score, the seven well-being items, the nine value-domain selections and the routed ratings Q11_3 to Q11_9, each in the coding the models use. `codebook_v53.csv` in this repository describes every column. The data file is not in this repository; how to obtain it is given in the Data Availability Statement of the article.

`00_derive.R` records how these variables were derived from the raw survey export (`data/data.xlsx`, sheets `data` and `Format`), which is not public. It is kept for that record and is not needed to run the analysis.

Every script is run from the repository root (`BASE <- "."`). The scripts create `outputs/` and write the derived dataset, the model caches and all tables and figures there or in the root.

## Software

R 4.5.1 with lavaan 0.6-21, poLCA 1.6.0.2, dplyr 1.1.4, readxl 1.4.5 (for `00_derive.R` only), car 3.1-3 and psych 2.6.3. Python 3.13 with matplotlib 3.10, Pillow 12.2 and python-docx 1.2. The full session record is in `sessionInfo.txt`.

## Run order

    Rscript 00_load_release.R
    Rscript 01_lca.R
    Rscript 02_structural.R 4 "Minimal,Relational-occupational,Family-health,Expansive"
    Rscript 03_value_change.R 4
    Rscript 04_breadth_distribution.R
    Rscript 04b_disturbance_fixing.R
    Rscript 06_breadth_group_ratios.R
    Rscript 07_reference_one_domain.R
    Rscript 08_disaster_modification.R
    Rscript 09_table_s12_panelA.R
    Rscript 10_sd_units_and_contrasts.R
    python 05_figures_python.py
    python 06_model_diagrams.py

- `00_load_release.R` builds `values_analysis.rds` from `data/data_release_v53.csv`, adding the variables derived from it (the number of endorsed domains, income in three bands for Table 1, and the mean ratings across endorsed domains), and runs the checks of the verification gate of `00_derive.R` that the data file allows. Any failed check stops the script before anything is written. `00_derive.R` builds the same dataset from the raw export, with the FAMCARE-10 score and the panel key in addition.
- `03_value_change.R` also tabulates the change ratings by FAMCARE-10 tertile when the FAMCARE-10 score is present. The data file does not carry it, so that tabulation is skipped; no table of the article uses it.
- `01_lca.R` fits the one- to six-class latent class models with 100 random starts each. The full enumeration takes more than an hour. The fitted models are saved to `lca_models_raw.rds` as soon as fitting ends, and a later run with the same settings reuses that file.
- `02_structural.R` needs the number of classes and their labels as arguments, and `03_value_change.R` needs the number of classes; neither chooses it. 02 fits an ordinal latent class model for a sensitivity check and caches it in `lca_ordinal_K4.rds`.
- `06_breadth_group_ratios.R` must run before `07_reference_one_domain.R`, because both write `outputs/table21_firstorder_breadth_tests.csv` and the article uses the version written by 07.
- `09_table_s12_panelA.R` evaluates lines 1 to 97 of `07_reference_one_domain.R` by line number. `10_sd_units_and_contrasts.R` evaluates the same lines and further blocks of 07, 02 and 08, located by their first and last lines. Editing those three scripts can therefore change what 09 and 10 compute.
- 09 and 10 fit no new model. They re-estimate the models defined in 02, 07 and 08 with the same syntax and stop unless the estimates reproduce the earlier outputs.
- Tables are rounded half up, as in the article: a value exactly half way goes up (63/112 = 56.25% is written 56.3), where base R's `round()` would go to the even digit (56.2). `00_load_release.R` and 02 to 10 source `round_half_up.R`, which makes `round()` and `signif()` in those scripts round half up; every `sprintf("%.Nf")` that writes a table cell is given a value already rounded that way. In 07 the `source()` call shares a line with `setwd(BASE)` so that the line numbers 09 and 10 rely on do not move. `00_derive.R` and `01_lca.R` do not source it; 01 writes its tables at full precision.

On one Windows machine with an ARM processor, `Rscript` returned exit code 127 at exit for most scripts, whether they had finished or stopped. There, the exit code cannot be used to tell success from failure. A script that stops prints `Error` and `Execution halted`; a script that finishes prints its own closing message.

## Where each table and figure comes from

All files are under `outputs/`. Figures are under `outputs/figures_plos/`, as TIFF files (LZW compression, RGB, 600 dpi, 7.5 inches wide), the format PLOS ONE requires.

| Article | Script | Output files |
|---|---|---|
| Table 1 | 02, and 07 for the rows on the number of endorsed domains | `table1_background.csv`, `table26_breadth_by_class.csv` |
| Table 2 | 02 | `table2_class_profile.csv` |
| Table 3 | 07 (coefficients and Wald tests), 10 (SD units) | `table22_breadth_ref1.csv`, `table22b_breadth_ref1_precise.csv`, `table23_ref1_wald.csv`, `table37_sd_units.csv`, `table38_sd_wb.csv` |
| Table 4 | 07, 10 | as Table 3, rows for the two family cancer groups |
| Table 5 | 08, 10 | `table31_disaster_group.csv`, `table31b_disaster_precise.csv`, `table32_disaster_wald.csv`, `table32b_disaster_wald_precise.csv`, `table37_sd_units.csv` |
| Fig 1 | none | Drawn in PowerPoint and exported as `Fig1.tif`; not produced by these scripts, and the PowerPoint file is not included |
| Fig 2 | 06_model_diagrams.py, `flow_diagram()` | `Fig2.tif`; the diagram carries no numbers |
| Fig 3 | 05_figures_python.py, `figure_endorsement()` | `Fig3.tif`, `Fig3_series.csv`, from `table2_class_profile.csv` and `table27_endorsement_rates.csv` |
| Fig 4 | 05_figures_python.py, `figure_coefficients()` | `Fig4.tif`, `Fig4_series.csv`, from `table37_sd_units.csv`, `table38_sd_wb.csv`, `table22_breadth_ref1.csv` and `table31_disaster_group.csv` |
| Table A in S1 Appendix | none | Questionnaire wording, compiled from the `Format` sheet of the raw export and the questionnaire draft, not by these scripts |
| Table B in S1 Appendix | 00 (Panel A), 02 (Panel B) | `gate_00_load_release.txt` or `gate_00_derive.txt`, `table15_reliability.csv` |
| Table C in S1 Appendix | 00 | the derivations in `00_derive.R` and `codebook_v53.csv`; `gate_00_derive.txt` (every check, from the raw export) or `gate_00_load_release.txt` (the checks the data file allows), `domain_codebook.csv` |
| Table D in S1 Appendix | 01 | `lca_fit_table.csv`, `lca_enumeration.md` |
| Table E in S1 Appendix | 02 (Panel A), 07 (Panel B), 10 (SD units) | `table13_class_nsup_overlap.csv`, `results_phase2.md` (R-squared, overlap coefficient), `table28_class_given_breadth.csv`, `table28b_class_given_breadth_tests.csv`, `table37_sd_units.csv` |
| Table F in S1 Appendix | 02, 10 (SD units) | `table10_breadth_matched.csv`, `table37_sd_units.csv` |
| Table G in S1 Appendix | 07, 10 (SD units) | `table25_firstorder_ref1_std.csv`, `table25b_firstorder_ref1_fit.csv`, `table21_firstorder_breadth_tests.csv`, `table41_s7_stdlv.csv` |
| Table H in S1 Appendix | several | A count of the tests reported in the article, compiled outside these scripts from `table9_test_inventory.csv` (02), `table21_firstorder_breadth_tests.csv` (07), `table23_ref1_wald.csv` and `table28b_class_given_breadth_tests.csv` (07), and `table32_disaster_wald.csv` and `table34_fourgroup_wald.csv` (08) |
| Table I in S1 Appendix | 01 | `lca_profiles.csv`, rows with K = 5 |
| Table J in S1 Appendix | 02 | `table11_nsup_ic.csv` |
| Table K in S1 Appendix | 03 | `phase3_by_domain.csv` |
| Table L in S1 Appendix | 09 (Panel A), 10 (Panel A, SD units), 07 (Panel B) | `table36_s12_panelA_breadth.csv`, `table37_sd_units.csv`, `table30_disturbance.csv` |
| Table M in S1 Appendix | 08 | `table33_fourgroup.csv`, `table34_fourgroup_wald.csv`, `table34b_fourgroup_wald_precise.csv`, `table35_cells.csv` |

`table37_sd_units.csv` has a `table` column that names the article table each row belongs to.

Outputs that do not appear in the article record checks made during the analysis: `table3_modelA.csv` to `table8_groupratio_by_predictor.csv`, `table12`, `table14` and the contrast and sensitivity files from 02; the other files from 03; `table16` from 04; `table17` from 04b; `table18` and `table20` from 06_breadth_group_ratios.R; `table24` and `table24b` from 07, which confirm that changing the reference category reparameterises the model without changing it; and `table39_contrasts.csv` and `table40_stdlv_delta.csv` from 10.

## SD units

Tables 3 to 5 and Tables E, F and L in S1 Appendix give each coefficient on the overall well-being score and also in SD units. SD units divide the coefficient, and its confidence limits where shown, by one number: the model-implied standard deviation of overall well-being in the primary model, `sqrt(lavInspect(fit, "cov.lv")["WB", "WB"])`. The same divisor is used for every model and every group. `table38_sd_wb.csv` lists the model-implied SD of overall well-being in each model and group for comparison. Table G in S1 Appendix is different: it reports each of the two well-being components in its own model-implied SD, taken from `standardizedSolution(type = "std.lv")`.
