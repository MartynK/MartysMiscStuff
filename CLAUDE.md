# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Purpose

This is a personal "miscellaneous" R package repository that serves as a bundled collection of standalone mini-projects, analysis scripts, and utilities. While structured as an R package (named `HUPHAR2024.PK.presentation` in DESCRIPTION, marked as "INCOMPLETE DEFUNCT"), it primarily functions as a working repository where each subfolder in `inst/` contains a distinct project or analysis.

**Core concept**: The `inst/` directory is the heart of this repository, containing independent projects that are kept together for organizational convenience rather than being true package components.

## Mini-Projects in inst/

### Shiny Applications

#### `BE_GMR_Bayes/`
**Purpose**: Bayesian bioequivalence analysis using GMR (Geometric Mean Ratio)
- Interactive Shiny apps for Bayesian BE assessment with different priors
- Uses `rstan` for Bayesian inference
- Precompiled Stan models stored in `comp_models/` (platform-specific: Windows vs Ubuntu)
- Three prior distributions: uniform, extreme (normal), wide_beta
- Evolution: `iter1.r` → `iter2_app.r` → `iter3_app.r` (most recent)
- **Running**: `setwd("inst/BE_GMR_Bayes"); shiny::runApp("iter3_app.r")`

#### `iscv_pooling/`
**Purpose**: Bayesian hierarchical pooling of intrasubject CVs from multiple studies
- Uses `cmdstanr` (not rstan)
- Beta(2,2) prior on CV, Inverse-Gamma likelihood
- Stan model: `stan/iscv_pool.stan`
- Structure: `server.R`, `ui.R`
- **Running**: `setwd("inst/iscv_pooling"); shiny::runApp()`

#### `cutoff_Shiny.app.R`
**Purpose**: Interactive visualization of acceptance/rejection regions
- Explore regions defined by `a*x^c + b*y^c > cutoff`
- Standalone single-file Shiny app

### Statistical Analysis Projects

#### `elephant_interim_improve/`
**Purpose**: Safety analysis for "ELEFANT" clinical trial
- RMarkdown report (`safety_v1.Rmd`) with extensive safety tables
- Author: Farkas Nelli
- Uses gtsummary, kableExtra for reporting
- Includes Excel data files

#### `vaers_query/`
**Purpose**: Pharmacovigilance disproportionality analysis using VAERS data
- Computes ROR (Reporting Odds Ratio) and PRR (Proportional Reporting Ratio)
- Queries VAERS Socrata API or works with local CSV files
- MedDRA Preferred Term matching
- Well-documented utility functions

#### `crossover_variabilities_plot/`
**Purpose**: Visualization options for 2×2 crossover trial data
- Creates paired dot-plots and alternative visualizations
- Demonstrates between-subject vs within-subject variability
- Uses patchwork for multi-panel figures

### Sample Size & Power Calculation Scripts

#### `dapparosu_sampsize.r`
Sample size for average bioequivalence (ABE) using `PowerTOST`
- Varying CV and PE (point estimate)
- Accounts for dropouts
- Visualizes sample size vs CV

#### `ci_length_ava_ddi_sampsize.r`
Confidence interval length exploration for 2×2 BE trials
- Shows relationship between n, CV, PE, and CI width
- NOTE: at PE=1, if CI length ≤ 30%, then 36n + dropouts adjustment

#### `ntg_sampsize_snippet.r`
Sample size AND cost calculations for NTG (nitroglycerin) study
- 2×2×4 crossover design
- Incorporates bioanalytical sample costs (€25/sample)
- Clinical costs (€11,000/subject)
- References actual study ISCVs (45.7% for study, 40.57% and 53.67% for products)

### Randomization Tools

#### `randomizer.r`
Stratified permuted block randomization using `randomizeR` package
- Multi-center (4 sites: A, B, C, D)
- Two severity strata (severe, light)
- 3-arm allocation (ratio 1:1:2)
- Block size = 4
- Includes simulation to assess balance

#### `caret_random.r`
Stratified randomization using `carat` package (StrPBR function)
- SONOMIND trial allocation example
- 4 sites with configurable probabilities
- Severity stratification (moderate 60%, severe 40%)
- Three treatment arms (High, Low, Sham)

### Pharmacokinetic Modeling

#### `nlmixr2_1stattempt.r`
PBPK (Physiologically-Based PK) modeling using `nlmixr2`
- Compartmental model with regional blood flows
- Partition coefficients (KbBR, KbMU, KbAD, KbBO, KbRB)
- Hepatic clearance (CLint)
- Uses xpose and xpose.nlmixr2 for diagnostics

#### `mesi_PK_curves.r` and `mesi_PK_curves_up..speccase.r`
PK curve visualization and analysis scripts

### Utility Scripts

#### `calculate_f2_different_points.r`
f2 similarity factor simulation and visualization
- Compares dissolution profiles
- Monte Carlo simulation (1M iterations)
- Generates 90% CI for f2 values

#### `net_satisfyer_score_eval.r`
Net promoter/satisfaction score calculation
- Computes net score (% top - % bottom)
- Standard error and 95% CI using delta method

#### `pics_stuff/merge_json.r`
**Purpose**: Merge Google Photos JSON metadata into JPEG EXIF data
- Reads Google Takeout JSON files
- Extracts timestamps, GPS coordinates, face regions
- Writes to JPEG EXIF using ExifTool
- Face region support using MWG-RS XMP standard

### Reference Materials

#### `lmer_cheat_sheet.rmd`
Comprehensive reference for `lme4` package
- Linear mixed-effects models
- Comparison with `nlme` capabilities
- Author: Clay Ford

#### `print_gtsummary.qmd`
Example of printing gtsummary tables to Word documents
- Uses `huxtable` for font size control

#### `report/`
Quarto report templates with modular structure
- `report.qmd`: Main report
- `child1_preamble.qmd`, `child2.qmd`: Child documents
- `custom_ref_doc.docx`: Word template for styling

### Shared Utilities

#### `functions/load_stuff.r`
Standard setup script for analysis projects. Sets:
- knitr options: no comments, no echo, no warnings, 300 DPI PNG output
- gtsummary themes: NEJM journal style, compact
- Scientific notation threshold (scipen = 1)
- Default digits = 3
- UTF-8 encoding

## Package Components

### Exported Functions (R/)
- `load_all_Rdata()`: Load all .Rdata/.rdata files from a directory
- `martys_format_p()`: Format p-values for presentation (handles <0.0001, rounds appropriately)
- `parallel_envir_setup.r`: Parallel computing environment setup
- `while_window_clicked_wrapper.r`: Window click wrapper utility

## Development Workflow

### Standard R Package Development
Per README.md, update in this order:
```r
# 1. Edit README.Rmd, then:
devtools::build_readme()

# 2. Update NEWS file (project name)

# 3. Update DESCRIPTION (dependencies, etc.)

# 4. After defining new functions under /R:
devtools::document()

# 5. Build pkgdown site:
devtools::build_site()

# 6. Install package:
devtools::install()

# 7. Check package:
devtools::check()
```

### Running Individual Projects
Most scripts in `inst/` are designed to be run directly:
```r
# For standalone scripts:
source("inst/dapparosu_sampsize.r")

# For Shiny apps, navigate to their directory first:
setwd("inst/BE_GMR_Bayes")
shiny::runApp("iter3_app.r")
```

### Testing
```r
devtools::test()
```
**Note**: Tests still reference "PROJECTNAME" placeholder - intentional for template workflow

## Key Dependencies by Domain

### Bayesian/Stan Modeling
- `rstan`: BE_GMR_Bayes apps
- `cmdstanr`: iscv_pooling app (requires CmdStan installation)

### Pharmacokinetics
- `nlmixr2`, `xpose`, `xpose.nlmixr2`: PK/PD modeling
- `PowerTOST`: Sample size for bioequivalence

### Randomization
- `randomizeR`: Permuted block randomization
- `carat`: Covariate-adaptive randomization

### Reporting & Tables
- `gtsummary`: Clinical tables (NEJM theme used throughout)
- `kableExtra`: Enhanced table formatting
- `huxtable`: Table printing to Word
- `rmdformats`: Pretty HTML outputs (downcute theme)

### Data Manipulation
- `dplyr`, `tidyr`: Tidyverse core
- `ggplot2`, `patchwork`: Visualization

### Other
- `here`: Path management
- `lubridate`: Date handling

## Important Notes

1. **Repository naming**: DESCRIPTION says "HUPHAR2024.PK.presentation" but this is intentionally generic. The repo is really "MartysMiscStuff" - a collection of projects.

2. **Platform-specific Stan models**: BE_GMR_Bayes uses different precompiled model files for Windows vs Ubuntu to avoid recompilation delays. The app detects OS via `Sys.info()[["sysname"]]`.

3. **Working directory assumptions**: Many `inst/` scripts assume:
   - Either run from their own directory, OR
   - Use `here::here()` for path resolution

4. **Standard styling**: Projects consistently use:
   - gtsummary with NEJM theme
   - 300 DPI PNG outputs
   - UTF-8 encoding
   - Compact, publication-ready formatting

5. **Git workflow**: Per README, use 'devel' branch for development work

6. **Incomplete package**: The package infrastructure exists mainly for organization; the real work is in `inst/` mini-projects that are largely independent.

## Common Commands

```r
# Load standard analysis environment
source(here::here("inst", "functions", "load_stuff.r"))

# Load all .Rdata files from a directory
load_all_Rdata(here::here("inst", "report", "backend"))

# Format p-value for publication
martys_format_p(0.00009)  # returns "<0.0001"
```
