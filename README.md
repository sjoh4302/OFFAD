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
- **Load previous study** - Select this option to view the clustering results from a past dataset and adjust selection of OFF periods. Users will be presented with a file navigation window with which to select a previously saved **OFFDATA** variable. To learn more of this choice, skip to section E

# A) Loading data
After choosing to begin a new study, users will be prompted to enter a number of details to set up the clustering process
- **Dataset name** - Enter a name for this study
- **VSspec pathin** - Select the matlab workspace containing the scored vigilance states for this dataset. Each state should be represented by its own array and named as following: 'w'=wake, 'nr'=NREM, 'r'=REM, 'mt'=movement
- **Epoch length** - length of scoring epochs in seconds
- **pNe pathin** - Select the matlab workspace containing your pNe signals. The worskpace should contain each channel as a seperate 1xN array in which N is the sample length. Channels should be named in the format 'Ch1' = Channel 1
- **pNe Fs** - Sampling rate of the pNe signal
- **pNe units** - Amplitude units of pNe data
- **LFPpathin** - (OPTIONAL) Select the matlab workspace containing your LFP signals. The worskpace should contain each channel as a seperate 1xN array in which N is the sample length. Channels should be named in the format 'Ch1' = Channel 1
- **pNe Fs** - (OPTIONAL) Sampling rate of the LFP signal
- **pNe units** - (OPTIONAL) Amplitude units of LFP data
- **Filter LFP** - If LFP used, select whether or not to bandpass filter the signal in the 0.5-100Hz range
- **Clustering evaluation method** - Choose which clustering evaluation method to use to select the optimum number of Guassian components (Maximum = 8 components). Daives-Bouldin and Calinsky-Harabasz methods are available. Optimal component number is computed from a random 0.1% sample of NREM data
- **Channels to remove** - Input the numbers of channels which should be ignored from the OFF-period detection pipeline, for example noisy channels. [1 2 3] = ignore channels 1, 2 and 3
-**Clustering variables** - Select the two variables with which to cluster to identify OFF-periods. Default is to use the time series of the smoothed pNe amplitude. Using this option, different smoothing strengths should be chosen for each variable. The second option is to use spectral power of the pNe signal as a variable (NOT RECOMMENDED)
-**Clustering sample size** - The size of the OFF-period component is computed from a subset of the total data consisting only of pNe data from NREM, then applied to all vigilance states. Choose the percentage of NREM pNe data to use for cluster size selection. A higher precentage will result in a longer clustering phase.

Once all fields have been correctly filled, select **Done** at the lower right corner of the GUI

# B) Preclustering
A new window will appear. This is the preclustering window which gives the user a indication of the OFF-period clustering on a subset of the data. A dropdown menu at the top shows the channel currently being analysed. Select which channel to view here. Two subplots will be displayed on the figure. The left subplot shows a density heatmap of the data for channel X in clustering variables 1 and 2. Lighter colours represent areas of high density. Channels with a good signal:noise ratio should show a clear cluster near the origin which represents the low amplitude OFF-periods. The right subplot shows the Gaussian Mixtue Model (GMM) generated from the data in the left plot with n components where n is the optimal number suggested by the clustering evaluation method. The component with the centroid with the lowest values in cluster variables 1 and 2 is shown in red an represents the OFF-period component. If the the red component on the right plot matches the dense OFF-period area of the left subplot, clustering should be succesful. If the clusters do not match, try re-clustering the data. Once you are happy, move to a new channel using the dropdown menu. Once you are happy, select **Done** to move on to the full clustering stage. The optimal component number computed for any of the channels in the preclustering stage will be carried forward. 

# C) Clustering
A warning window will appear to inform the user that clustering is taking place. Each channel is clustered seperately. The number of GMM components (n) will be the number identified in the preclustering stage. Any channels not viewed in th preclustering stage will have their optimal GMM component number determined automatically. A GMM model with n components will then computed from a random sub-sample of the NREM data, the size of which was chosen by the user. Next, the full dataset will be assigned into clusters based on their component posterior probability (i.e likelhood of belonging to each component). Data points belonging to the component with the lowest mean in clustering variables 1 and 2 are classified as OFF-period points. All other data points belonging to all other components are classified as ON-period points. 

# D) Summary statistics
