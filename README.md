# CIBI-biodiversity
Scripts and workflows used for data processing, analyses, and visualization for the California Insect Biodiversity Initiative (CIBI) field collections manuscript titled: DNA barcoding of more than one million insect specimens indicates that at least one third of California’s insect biodiversity remains undiscovered.

Note that the input data file `250626_BOLD_data_clean.xlsx` is too large to be hosted on Github, and can be obtained from the Supplementary Material from the manuscript. 


### Read me for BIN Distribution Map Pipeline ###

### Step 1: Download data files

## Step 1a: Download BOLD Excel Data Spreadsheets for the project from boldsystems.org
- Data files are split into subprojects of ~50k specimens to download smaller files
- Select 'Progress Report', 'Taxonomy', and 'Collection Data' as sheets to include, Multi-Page format

## Step 1b: Download LandFire datasets
- For California, download CONUS Existing Vegetation Cover for the most recent year (2023 for present study)
- https://landfire.gov/data/FullExtentDownloads

## Step 1c: Download Ecoregion dataset
- Download Level III and Level IV Shapefiles
- https://www.epa.gov/eco-research/ecoregion-download-files-state-region-9#pane-04

### Step 2: Run data cleaning script (can do on local computer)

- california_grid_landfire_data_prep.R
- this requires manual fixes about halfway through
- this script generates the following files: 
- filtered_data.txt: combined all of the BOLD spreadsheets and filtered taxonomically and geographically
- filled_data.xlsx: datasheet that includes easy taxonomic corrections, meant for manual correction for custom problems; save as filled_data_clean.xlsx once finsihed
- CA_landfire_raster.tif: save time by re-uploading this file rather than generating it again for the landfire raster
- CA_lf_p.tif: the final landfire + ecoregion raster that will be used in downstream analysis
- CA_levels.rds: this file contains names and numbers associated with landfire vegetation types
- BOLD_data_lf.rds: this is the final specimen data file with geospatial data added

### Step 3: Transfer files to cluster

- Copy California_landfire_grid_template folder into working directory
- Copy these files into new folder's data subfolder: BOLD_data_lf.rds, CA_lf_p.tif, CA_levels.rds
- Copy data folder into 'big' and 'small' analysis folders

### Step 4: Generate BIN maps
- Go to the California_grid_landfire_big folder
- Run 2a_California_grid_landfire_cluster_big.R 
- Go to the California_grid_landfire_small folder
- Run 2b_California_grid_landfire_cluster_small.R
- Outputs: folders in each top folder (big and small) 'ecoregion_maps' and 'vegtype_ecoregion_maps' full of maps for downstream analysis

### Step 5: Run scripts to resize and combine maps
- Run scripts 3–7 to resize maps and combine maps for groups by taxon

