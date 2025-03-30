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
   Scripts/01_read_data.R
   Scripts/02_search_aglycon.R
   Scripts/03_match_msms.R
   # optional
   Scripts/03b_search_msms_negativemode.R 
   ```
## Output
- Tables: "output/"
- MS/MS matched results: "output/db/"

## Reproducibility
- R version: See "session_info.txt"
- Dependencies: See "session_info.txt"; Scripts/00_setup.R
  