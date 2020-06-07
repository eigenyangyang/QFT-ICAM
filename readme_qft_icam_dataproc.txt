Python codes for processing data measured by QFT-ICAM (Quantitative Filter 
Technique - Integrating Cavity Absorption Meter) system. 
Script name: qft_icam_postprocess.py

------------
Preparations
------------
Before using the scripts, please:
1 install:
- Python 3.6x or higher version    
- Essential Python packages: numpy, scipy, pandas, matplotlib.

2 correct labelling in .raw files according to the protocol.
1) Use uniform labeling method for all samples. For example, the 2 measurements of the sample "UW5" should be labeled as "UW5#1" and "UW5#2" (here separator is '#' and should be used throughout all samples) or simply both as "UW5". 
2) Labels for blank filters measurements for contamination control of each box of GF/F filters (they are NOT reference measurements!) need to contain keywords "batch", "MQ" or "SW" (case insensitive), e.g."*batch_MQ#11". Normal sample labels should not contain any aforementioned keywords.

3 modify the path of the working directory according to your situation in the main program code block (Lines below "if name == 'main':"). 

4 prepare the following files (see examples and data processing procedure below) and put them in the same directory as the one containing .qft.raw files.
- config_PostProc.txt
- rmSpectraIndex.txt
- qft_icam_matched_labels_filtration_volumn.txt
- labels_latlon_datetime_qft.txt


----------------
Data Processing
----------------

-------------------------------------------------------------------
Step 1: Calculate optical density spectra (OD) from .qft.raw files.
-------------------------------------------------------------------
a) Change the working directory as the one storing .qft.raw data files 
by modifying the variable "wd" (Line ~887).
b) Open the script, comment out Step 2-5 with "#", and then save file.
c) In the terminal, run the script by:
############################################
python3 qft_icam_postprocess.py
############################################

Step 1 calculates sample OD by correcting signals from stray light, dark current, reference filter, diffuse light field inside the ICAM, and chlorophyll fluorescence (not considered for bleached samples). For each measurement, OD is calculated twice, with reference before and after sample, respectively.

Output: Level 0 OD data (.l0a files).

-------------------------------------------------------------------
Step 2: Plot each sample OD in a .l0a file.
-------------------------------------------------------------------
a) Open the script and comment out Step1, 3-5.
b) Run the script.

Output:
Folders named as the base name of the .l0a file in the working direcoty 
containing plots from individual sample measurements. The legend of each 
absorption spectrum in every plot gives the row index of each data point 
when this .l0a data file is imported into a pandas dataframe via "pd.read_csv". 
These row indices can be used as input of Step 3 to remove suspicious data.

--------------------------------------------------------------------------
Step3.1: Manuelly remove suspicious measurements in .l0a data indexed by 
"row_index" (indicated in the plot legend, see Step 2) and take the median 
values and standard deviation of repeatedly measured data.

Step3.2: Merge all '*_OD_median.txt' and '*_OD_sd.txt' files, seperates total particles 
(sample), non-algal particles (bleached sample) and contamination control blank 
filter measurements, and merges data from each catalog, respectively.
--------------------------------------------------------------------------
a) Open the script and comment out Step1-2, 4-5.
b) Prepare your own "rmSpectraIndex.txt" file (tab delimited) based on the plots generated from Step2 in the working directory.
c) Run the script.


Output of Step3.1:
1) Tab delimited files in the direcoty "ODfiles" named as "*_OD_median.txt" and "*_OD_sd.txt", with "*" standing for the base name of the .l0a file in the working direcoty,respectively.
2) Plots of the spectra in the "*_OD_median.txt" and "*_OD_sd.txt" named as "*_OD_median.png" and "*_OD_sd.png", respectively, in the plotting directory (see Step2). 


Output of Step3.2:
1) Eight tab delimited files in the directory "ODfiles" named as 
    "qft_icam_merged_median_OD_all.txt" (OD spectra from all measurements), 
    "qft_icam_merged_sd_OD_all.txt" (standard deviation of OD spectra from all measurements), 
    "qft_icam_merged_median_OD_totparticle.txt" (OD spectra from samples), 
    "qft_icam_merged_sd_OD_totparticle.txt" (standard deviation of OD spectra from samples), 
    "qft_icam_merged_median_OD_NAP.txt" (OD spectra from bleached samples),
    "qft_icam_merged_sd_OD_NAP.txt" (standard deviation of OD spectra from bleached samples),
    "qft_icam_merged_median_OD_batch.txt" (OD spectra from blank filters), and
    "qft_icam_merged_sd_OD_batch.txt" (standard deviation of OD spectra from blank filters), respectively.
2) Four plots named as 
    "qft_icam_merged_median_OD_totparticle.png" (OD spectra from samples).
    "qft_icam_merged_sd_OD_totparticle.png" (standard deviation of OD spectra from samples).
    "qft_icam_merged_median_OD_NAP.png" (OD spectra from bleached samples), 
    "qft_icam_merged_sd_OD_NAP.png" (standard deviation of OD spectra from bleached samples), respectively.
3) A file named as "qft_icam_merged_labels.txt" containing all the data lables.
This can be used to construct the files "qft_icam_matched_labels_filtration_volumn.txt" and "labels_latlon_datetime_qft.txt". Note that the sample and bleached sample labels are changed to "sam_" or "BL_" + original labels.

--------------------------------------------------------------------------
Step 4: Calculate absorption coefficients of total particulate 
matters, non-algal particles and phytoplankton.
--------------------------------------------------------------------------
a) Open the script and comment out Step1-3, 5.
b) Prepare your own "qft_icam_matched_labels_filtration_volumn.txt" file (tab delimited) based on the file 'qft_icam_merged_labels.txt' generated from Step3.2 in the working directory.
c) Run the script.

Output:
1) Absorption coefficients of total particulate matters and NAP, and their standard deviations saved in the working directory as 
"qft_icam_merged_median_abs_totparticle.txt",
"qft_icam_merged_median_abs_NAP.txt",
"qft_icam_merged_sd_abs_totparticle.txt", 
"qft_icam_merged_sd_abs_NAP.txt", respectively. 
The first 2 columns are sample labels and filtration volumn. 
2) Adjusted NAP absorption coefficient data saved in the working directory as "qft_icam_merged_median_abs_NAP_adjust.txt". "Adjusted" means bringing 
NAP absorption back to ap in the NIR so that aph in NIR is 0 (offset the 
average values in 720-750nm) (Neeley et al. 2018).
3) Adjusted NAP absorption coefficient data plotted in the working directory
as "qft_icam_merged_median_abs_NAP_adjust.png".
4) Phytoplankton absorption coefficient data calculated using adjusted NAP absorption and saved in the working directory as "qft_icam_merged_median_abs_phytoplankton.txt". 
5) Phytoplankton absorption coefficient data plotted in the working directory as "qft_icam_merged_median_abs_phytoplankton.png".
6) Each ap, a_NAP and aph spectrum plotted in the directories of "abs_totparticle", "abs_NAP", "abs_phytoplankton", respectively. 

--------------------------------------------------------------------------
Step 5: Prepare data for upload to Pangaea.
--------------------------------------------------------------------------
a) Open the script and comment out Step1-4.
b) Prepare your own "labels_latlon_datetime_qft.txt" file (tab delimited) in the working directory.
c) Run the script.

Output:
QFT-ICAM data suitable to be uploaded to Pangaea.
1) Five tab delimited absorption (and standard deviation) files in the 
working directory named as
"qft_icam_merged_median_abs_totparticle_pangaea.txt",
"qft_icam_merged_sd_abs_totparticle_pangaea.txt",
"qft_icam_merged_median_abs_NAP_adjust_pangaea.txt",
"qft_icam_merged_sd_abs_NAP_pangaea.txt",
"qft_icam_merged_median_abs_phytoplankton_pangaea.txt".
2) Four tab delimited OD (and standard deviation) files in the directory 
"ODfiles" named as
"qft_icam_merged_median_OD_totparticle_pangaea.txt",
"qft_icam_merged_sd_OD_totparticle_pangaea.txt",
"qft_icam_merged_median_OD_NAP_pangaea.txt",
"qft_icam_merged_sd_OD_NAP_pangaea.txt".

----------
References
----------
1) Röttgers et al. (2016). Quantitative filter technique measurements of spectral light absorption by aquatic particles using a portable integrating cavity absorption meter (QFT-ICAM). Optics express, 24(2), A1-A20.
2) IOCCG Protocol Series (2018). Inherent Optical Property Measurements and Protocols:Absorption Coefficient, Neeley, A. R. and Mannino, A. (eds.), IOCCG Ocean Optics and Biogeochemistry Protocols for Satellite Ocean Colour Sensor Validation, Volume 1.0, IOCCG, Dartmouth, NS, Canada.
