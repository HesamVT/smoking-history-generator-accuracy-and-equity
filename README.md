# Smoking History Generator Accuracy and Equity Research Archive

Hesam Mahmoudi  
Lung Cancer Policy Team (LCP)  
MGB Center for Health Technology Assessment (CHTA)  
Massachusetts General Hospital (MGH)  
Harvard Medical School  
Prepared: 2026-05-04  
Updated: 2026-09-24

This folder contains the analytic code and literature inventory prepared for public sharing alongside the manuscript appendix. The archived code evaluates the accuracy of Smoking History Generator-simulated smoking histories, their implications for lung cancer screening eligibility, and related equity concerns across SCCS, MECS, and BRFSS. The literature inventory documents peer-reviewed publications that used the National Cancer Institute Smoking History Generator (SHG), its outputs, or an explicit adaptation or derivative of its methods.

Included files:

- `SHG SCCS.Rmd`
- `SHG MECS.Rmd`
- `SHG BRFSS.Rmd`
- `shg_aggr_results.py`
- `SHG_literature_inventory_2026.xlsx`
- `SHG_literature_inventory_2026.csv`

## Scope

The three R Markdown files contain the cohort-specific analysis code for:

- SCCS
- MECS
- BRFSS

The Python file contains the aggregate post-processing code that combines completed cohort-level exports into manuscript-ready summary tables and figures aligned with the manuscript's cross-cohort accuracy and screening-eligibility analyses.

## SHG Literature Inventory

The literature inventory contains 143 unique peer-reviewed publications identified through September 2026. Publications are ordered by year and include the following fields:

- row number
- author(s)
- publication year
- full title
- journal
- DOI, when available
- SHG engagement classification
- paper-specific justification for inclusion

The inventory covers complementary publication pools, including citations to foundational SHG methods and input-development papers; applications from CISNET lung cancer modeling groups and established CISNET model families; publications by investigators associated with these models; subgroup, state, and regional extensions; international adaptations; tobacco-control and smoking-and-vaping models derived from SHG inputs or methods; and applications outside lung cancer. Publications that only cited or discussed SHG without using its methods, inputs, outputs, or documented model lineage were excluded.

The `SHG engagement` field describes the primary form of connection to SHG. It distinguishes direct use in an analysis; use of SHG-derived histories, transition rates, intensity distributions, mortality inputs, or calibration inputs; use through a documented CISNET model lineage; development or synthesis of the SHG framework and inputs; geographic or population-specific adaptations; explicit reimplementations; and derivative population or tobacco-policy models. This classification is descriptive and is not a rating of study quality.

`SHG_literature_inventory_2026.xlsx` is the formatted version of the inventory. `SHG_literature_inventory_2026.csv` contains the same publication-level fields in a machine-readable format that can be viewed and searched directly on GitHub.

## Aggregate Script

The aggregate workflow is implemented in `shg_aggr_results.py`.

It reads precomputed cohort-level outputs and renders:

- aggregate median tables
- aggregate threshold tables using ever-smoker denominators
- aggregate threshold tables using total-population denominators
- aggregate probability tables for the `P` and `2P - 1` summaries
- the combined eligibility Figure 1

The aggregate script is a post-processing step. It does not recompute the cohort-level analyses from raw data.

## Materials Not Included

The following materials are required to run the full workflows but are not included in this shared code archive:

- SHG simulator source, binaries, and auxiliary input files
- restricted cohort data files
- intermediate cohort-generated CSV, HTML, image, and analysis-workbook outputs

## Cohort-Specific Requirements

Each cohort R Markdown file depends on:

- the relevant cohort smoking-history dataset
- the SHG simulator files under `SHG6.3.4`
- a Windows executable for SHG runs where the code calls `lbc_smokehist_win.exe`
- the R packages loaded in each file header

For SCCS specifically, the code reads the restricted cohort data file:

- `SCCS_with_smoking_info_Lung_Breast_Colon.dta`

For MECS specifically, the code reads the restricted cohort data file:

- `MEC_Small.dta`

For BRFSS specifically, the code reads the restricted cohort data file:

- `BRFSS_Data_Small.dta`

## Aggregate Script Requirements

The aggregate Python script depends on completed cohort-level exports produced upstream by the cohort workflows, including:

- median table CSVs
- threshold table CSVs
- probability table CSVs
- eligibility heatmap HTML tables

These inputs are not included here because they derive from restricted cohort data and non-shared SHG resources.

## Notes

- The R Markdown files include their own inline narrative and run instructions.
- The aggregate script expects completed cohort-generated CSV and HTML inputs to exist before it is run.
- This archive is intended to document the analytic code structure used in the manuscript, not to serve as a fully self-contained reproducibility package.
