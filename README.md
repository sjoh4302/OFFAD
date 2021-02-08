# OFFAD (Off period automated detection)

# What is OFFAD?
OFFAD is an open source signal processing GUI for detecting OFF periods: periods of generalized neuronal silence associated with NREM slow waves. The pipeline consists of a Guassian clustering step followed by a visual and statistical presentation of the resulting OFF periods to allow further curation by the user.  

# Installing
You can clone this repository or download a ZIP file by clicking the **Code** button.

# Starting OFFAD
1. Open MATLAB 
2. Navigate to the **OffPeriodDetection** folder containing the OFFAD scripts
3. Run the following in MATLAB
```
OFF_AD
```
4. Alternatively, open the **OFF_AD.m** script and click **Run**

# Using OFFAD
Upon opening OFFAD, users will be presented with the following choices
- **Open new study** - Select this option to begin clustering a new dataset
- **Load previous study** - Select this option to view the clustering results from a past dataset and adjust selection of OFF periods. Users will be presented with a file navigation window with which to select a previously saved **OFFDATA** variable. To learn more of this choice, skip to section D

# A) Loading data
After choosing to begin a new study, users will be prompted to enter a number of details to set up the clustering process
- **Dataset name** - 

