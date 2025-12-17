# River Otter Species Distribution Modeling in Texas  
**GBIF → Maxnet → Scenario Projections**

**Author:** ANONYMOUS FOR REVIEW  
**Affiliation:** ANONYMOUS FOR REVIEW 

**Correspondence:** ANONYMOUS FOR REVIEW  

---

## Overview

This repository contains the data, code, and derived outputs used to model the
distribution of North American river otters (*Lontra canadensis*) across Texas,
USA. Species distribution models (SDMs) were developed using occurrence records
from the Global Biodiversity Information Facility (GBIF), hydrological predictors
derived from the National Hydrography Dataset (NHD), and the Maxent algorithm
implemented via the `maxnet` R package.

The repository is structured so that a reviewer or reader can:
1. Reproduce all analytical steps reported in the manuscript,
2. Regenerate intermediate products if desired, and
3. Inspect all manuscript-reported tables and figures.

Connectivity and cartographic post-processing were conducted in ArcGIS Pro, as
described in the manuscript.

---
---

## Reproducibility Instructions

This repository is designed to support full computational reproducibility.

1. Clone this repository (or download as a ZIP):
   ANONYMOUS FOR REVIEW

2. Download large raster datasets from Zenodo (DOI: XXXXX).

3. Place the downloaded files into the following directories:

   - `00_raw/nhd/`
     - `natwatfocal.tif`
     - `humwatfocal.tif`
     - `allwatfocal.tif`

   - (Optional, if not regenerating predictors)
     - `01_intermediate/`
       - `dens_all_km2.tif`
       - `dens_nat_km2.tif`
       - `dens_hum_km2.tif`
       - `dist_major_rivers_km.tif`

4. Open `RiverOtterSDM.Rproj` in RStudio.

5. Run `SDM_Pipeline_EndToEnd.R`.

Set `overwrite = TRUE` to regenerate all intermediate products and predictions,
or `overwrite = FALSE` to reproduce manuscript-reported results from existing
outputs.

## Data Sources

### Occurrence Data
- **GBIF**: Global Biodiversity Information Facility  
  - Species: *Lontra canadensis*  
  - Geographic scope: Texas, USA  
  - Temporal scope: 1900–2025  
  - GBIF records were filtered for coordinate validity, temporal completeness,
    and duplicate locations prior to analysis.

### Environmental Predictors
- **National Hydrography Dataset (NHD)**  
  - Natural water features  
  - Anthropogenic (human-modified) water features  
- **Major Rivers (Texas)**  
- **Hydrologic Unit Code (HUC8) watersheds**

All spatial processing was conducted using an equal-area projection
(NAD83 / CONUS Albers, EPSG:5070) where appropriate.

---

## Repository Structure

The folder structure mirrors the analytical pipeline and is referenced directly
by the R scripts. Folder names should not be changed.

# River Otter Species Distribution Modeling in Texas  
**GBIF → Maxnet → Scenario Projections**

**Author:** ANONYMOUS FOR REVIEW  
**Affiliation:** ANONYMOUS FOR REVIEW

**Correspondence:** ANONYMOUS FOR REVIEW 

---

## Overview

This repository contains the data, code, and derived outputs used to model the
distribution of North American river otters (*Lontra canadensis*) across Texas,
USA. Species distribution models (SDMs) were developed using occurrence records
from the Global Biodiversity Information Facility (GBIF), hydrological predictors
derived from the National Hydrography Dataset (NHD), and the Maxent algorithm
implemented via the `maxnet` R package.

The repository is structured so that a reviewer or reader can:
1. Reproduce all analytical steps reported in the manuscript,
2. Regenerate intermediate products if desired, and
3. Inspect all manuscript-reported tables and figures.

Connectivity and cartographic post-processing were conducted in ArcGIS Pro, as
described in the manuscript.

---

## Data Sources

### Occurrence Data
- **GBIF**: Global Biodiversity Information Facility  
  - Species: *Lontra canadensis*  
  - Geographic scope: Texas, USA  
  - Temporal scope: 1900–2025  
  - GBIF records were filtered for coordinate validity, temporal completeness,
    and duplicate locations prior to analysis.

### Environmental Predictors
- **National Hydrography Dataset (NHD)**  
  - Natural water features  
  - Anthropogenic (human-modified) water features  
- **Major Rivers (Texas)**  
- **Hydrologic Unit Code (HUC8) watersheds**

All spatial processing was conducted using an equal-area projection
(NAD83 / CONUS Albers, EPSG:5070) where appropriate.

---

## Repository Structure

The folder structure mirrors the analytical pipeline and is referenced directly
by the R scripts. Folder names should not be changed.

RiverOtterSDM_Paper/
├── 00_raw/ # Raw, unmodified input data
├── 01_intermediate/ # Derived predictor rasters
├── 02_occ/ # Cleaned occurrences and extracted covariates
├── 03_models/ # Model objects, tuning results, and coefficients
├── 04_pred/ # Scenario prediction rasters and area summaries
├── 05_figs/ # Manuscript figures
├── 06_connectivity/ # GIS-based connectivity inputs and outputs
├── 07_results/ # Final tables used in the manuscript
├── SDM_Pipeline_EndToEnd.R
├── RiverOtterSDM.Rproj
└── README.md

## Large Raster Files and Storage Policy

Due to GitHub file size limits, large raster datasets required to run the full SDM pipeline are archived on Zenodo. These include both raw hydrological focal rasters and derived intermediate predictor surfaces.

The following folders are included in this repository as placeholders but do not contain raster files:

00_raw/nhd/ — raw hydrological focal rasters
(natwatfocal.tif, humwatfocal.tif, allwatfocal.tif)

01_intermediate/ — derived predictor rasters
(water density and distance-to-river surfaces)

These rasters can be obtained from the Zenodo archive associated with this project:

Zenodo archive: [DOI to be added upon publication]

To reproduce the analysis without regenerating these rasters:

Download the Zenodo archive

Place the raw focal rasters into 00_raw/nhd/

Place the derived predictor rasters into 01_intermediate/

Run SDM_Pipeline_EndToEnd.R with overwrite = FALSE

Alternatively, setting overwrite = TRUE will regenerate all intermediate products from raw inputs.

---

## Folder Descriptions

### `00_raw/` — Raw Inputs (Do Not Modify)
Contains all original datasets used in the analysis, including:
- GBIF occurrence export
- NHD water feature layers
- Major rivers shapefile
- State boundary and watershed layers
- Lookup table used to classify NHD feature codes as natural or anthropogenic

These files are preserved in their original form for transparency.

---

### `01_intermediate/` — Predictor Rasters
Derived rasters created from raw hydrological layers, including:
- Natural water density (km²)
- Anthropogenic water density (km²)
- Combined water density (km²)
- Distance to major rivers (km)

These rasters are inputs to the SDM and can be regenerated by rerunning the pipeline.

---

### `02_occ/` — Occurrence Processing
Contains:
- Cleaned and filtered GBIF occurrence records
- Spatially thinned presence locations
- Extracted environmental covariates for presence and background samples

These tables define the final modeling dataset.

---

### `03_models/` — Model Outputs
Includes:
- K-fold cross-validation tuning results
- Fold-level performance metrics
- Final fitted Maxent model (`.rds`)
- Omit-10% training threshold
- Raw feature coefficients
- Predictor importance summaries (Σ|β|)

These files support model evaluation and interpretation.

---

### `04_pred/` — Scenario Predictions
Raster outputs used for mapping and spatial analysis, including:
- Current hydrology suitability (cloglog)
- No-human-water suitability (cloglog)
- Binary suitability rasters (thresholded)
- Difference rasters (current − no-human)
- Suitable area summaries (km²)

These rasters were used for GIS-based analyses and figure production.

---

### `05_figs/` — Figures
Final figures included in the manuscript, such as:
- Predictor importance plots
- Response curves
- Distribution and scenario maps
- Cost and accessibility summary plots

Intermediate or draft figures have been removed.

---

### `06_connectivity/` — Connectivity Analysis
Spatial layers and tables used for least-cost and patch-connectivity analyses,
including:
- HUC8 centroids (all and occupied)
- Watershed layers with otter occupancy
- Patch connectivity summary tables

Connectivity modeling was conducted in ArcGIS Pro.

---

### `07_results/` — Manuscript Tables and Logs
Contains final tables reported in the manuscript, including:
- County-level occurrence summaries
- First-detection period summaries
- Habitat suitability area comparisons
- Connectivity comparisons
- CostAccessibility.csv
Watershed-level accumulated cost distance metrics exported from ArcGIS Pro.
Columns include natural-only and contemporary accumulated cost distances.

- CostAccessibility_summary_stats.csv
Descriptive statistics (mean, SD, median, IQR, quantiles) for natural,
contemporary, and Δcost metrics.

- CostAccessibility_correlations.csv
Pearson and Spearman correlations among cost metrics.

- CostAccessibility_tests.txt
Paired statistical comparisons between natural and contemporary cost
distances, including effect sizes.

Also includes:

- `pipeline_log.txt`: execution log from the R pipeline
- `sessionInfo_full_pipeline.txt`: R session information for reproducibility

---

## Reproducing the Analysis

The full analysis can be reproduced by running:

SDM_Pipeline_EndToEnd.R

from the project root directory. The script is designed to run end-to-end and
prints key values used in the manuscript (e.g., sample sizes, thresholds,
predictor importance) to the console and log files.

An `overwrite` flag at the top of the script controls whether existing outputs
are reused or regenerated.

---

## Software Requirements

- R (≥ 4.2)
- Packages: `terra`, `sf`, `maxnet`, `dplyr`, `ggplot2`, `tigris`, `tidyr`,
  `readr`, `lubridate`, and dependencies listed in `sessionInfo_full_pipeline.txt`
- ArcGIS Pro (for connectivity and cartographic post-processing only)

---

## Citation

If you use these data or code, please cite the associated manuscript and the
original data sources, including GBIF and the National Hydrography Dataset.

---

## License

This repository is provided for academic and non-commercial research use.
Please cite appropriately and do not redistribute raw third-party datasets
without adhering to their original licenses.