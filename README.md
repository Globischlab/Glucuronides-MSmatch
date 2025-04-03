# Glucuronides_MSmatch
## PhaseII-MSmatch
This project is to build a streamline to identify the glucuronides and its aglycon after the enzymatic treatment, then perform the MS1 and MS2 library search for annotation.

## Usage
1. Clone this repo:
   ```bash
   git clone https://github.com/fyang22/Glucuronides-MSmatch.git
   ```
2. Install R dependencies:
   ```bash
   Scripts/00_setup.R
   ```
3. Run scripts sequentially:
   ```bash
   # read mzml data from data/raw
   Scripts/01_read_data.R
   # 1. compare the precursor masses in control samples and precursor masses after Enzymatic treatment: mz_aglycon = mz_control - 176.0321
   # 2. LC-MS precursor mass search with hmdb and pubchem substract glucuronides
   Scripts/02_search_aglycon.R
   # pre-select the tentative annotated mz and perform ms/ms screening
   Scripts/03_match_msms.R
   # Negative ionization mode only optional, search the characteristic fragmentations of glucuronides
   Scripts/03b_search_msms_negativemode.R
   # ms/ms screening with aglycon
   Scripts/04_match_msms_aglycon.R
   ```

## data
- databases: compounds list from hmdb or pubchem
- negative_ion_mode: raw data for negative ionization mode
- positive_ion_mode: raw data for positive ionization mode
  
## Output
- Tables: "output/"
  - spectra data
  - aglycons search in ea data
  - glucuronides and aglycon search in hmdb and pubchem
    
- MS/MS matched results: "output/db/"

## Reproducibility
- R version: See "session_info.txt"
- Dependencies: See "session_info.txt"; Scripts/00_setup.R
  