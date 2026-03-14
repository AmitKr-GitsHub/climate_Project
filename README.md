# Geospatial Climate Data and Economic Outcomes (India)

This repository contains assignment deliverables (excluding the in-class presentation):

- `code/climate_econ_assignment.R`: End-to-end R pipeline for ERA5 precipitation processing,
  climate shock construction, panel merge, and econometric estimation.
- `writeup/report.tex`: Concise LaTeX write-up of methodology and empirical strategy.
- `data/clean/economic_outcome_monthly.csv`: Small example cleaned economic file
  (schema template for merge reproducibility).

## How to run

1. Install R packages used in the script (`data.table`, `sf`, `terra`, `exactextractr`, `fixest`, etc.).
2. Place ERA5 netCDF file at:
   - `data/raw/era5/era5_tp_india_daily.nc`
3. Run:
   - `Rscript code/climate_econ_assignment.R`

## Notes

- The script includes an optional CDS API helper function for ERA5 download.
- Replace `data/clean/economic_outcome_monthly.csv` with your real outcome dataset,
  keeping columns: `unit_id`, `ym`, `y_outcome`.
