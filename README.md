# Smoking History Generation Code Archive

Hesam Mahmoudi  
Lung Cancer Policy Team (LCP)  
MGB Center for Health Technology Assessment (CHTA)  
Massachusetts General Hospital (MGH)  
Harvard Medical School  
Prepared: 2026-05-04

This folder contains the code files prepared for public sharing alongside the manuscript appendix.

Included files:

- `SHG SCCS.Rmd`
- `SHG MECS.Rmd`
- `SHG BRFSS.Rmd`
- `shg_aggr_results.py`

## Scope

The three R Markdown files contain the cohort-specific analysis code for:

- SCCS
- MECS
- BRFSS

The Python file contains the aggregate post-processing code that combines completed cohort-level exports into manuscript-ready summary tables and figures.

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
- intermediate cohort-generated CSV, HTML, image, and workbook outputs

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
