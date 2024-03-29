# Meerkat vocal interactions code

This repository contains the code to carry out analyses of meerkat vocal interactions,
described in the paper **"Mapping vocal interactions in space and time differentiates signal broadcast vs. signal exchange in meerkat groups"**.

The figures produced by running the code here include:

* Call category bigrams (Figure 1)

* Self-reply vs caller exchange (Figure 2)

* Tree plot (Figure 3)

* Call response dynamics curves for close calls and short notes as a function of distance between caller and responder (Figure 4a)

* Call response dynamics curves for close calls and short notes as a function of age (Figure 4b)

* Clustering of close calls and short notes over different spatial and temporal scales (Figure 5)

* Mean call rate vs speed and number of nearby neighbors (Figure 6)

## Bigrams and caller exchange (Figures 1, 2, and 3)

This section explains how to run the analyses to generate Figures 1 - 3 in the paper.

### Requirements - libraries

See the libraries required at the top of each script.

### Requirements - data

The code requires the following data files to run:

* `03_tree_data.RData`

* `04_meerkat_bigrams.RData`

* `05_self_nonself_reply_proportions.RData`

### Scripts

The scripts to produce the first 3 figures in the manuscript are the following:

* `03_min_tree_script.R`: produces Figure 3

* `04_meerkat_bigrams.R`: produces Figure 1

* `05_self_nonself_reply_proportions.R`: produces Figure 2

### How to run

Each script is associated with a data file of the same name. Run the scripts to produce the plots.

## Call-response and call clustering (Figures 4 - 5)

This section explains how to run the call-response dynamics analysis and the call clustering analysis to generate Figures 4 and 5 in the paper.

### Requirements - libraries

The code was developed and tested in R version 4.1.2 on Mac OS 13.2.1, Platform x86_64-appledarwin17.0 (64-bit). 
Other versions of R and operating systems are likely to work as well.
The code requires the following R packages to be installed via `install.packages('packageName')` prior to use:

* `viridis`: used for plotting

* `lubridate`: used for datetime handling

* `stringr`: used for string handling

* `gplots`: used for plotting

* `fields`: used for plotting


### Requirements - data

The code requires the following raw data files to run:

* `all_calls_sync_resolved_2023-09-10_cc_sn_filt_with_oor_anon.csv`: contains a table of all labeled calls

* `HM2017_COORDINATES_all_sessions_anon.RData`: contains GPS data (UTM) for the group HM2017 as well as timestamps

* `HM2017_DATAPRESENCE_all_sessions_anon.RData`: contains information on whether GPS and audio recordings were on at every second of the dataset for group HM2017

* `HM2019_COORDINATES_all_sessions_anon.RData`, `HM2019_DATAPRESENCE_all_sessions_anon.RData`, `L2019_COORDINATES_all_sessions_anon.RData`, `L2019_DATAPRESENCE_all_sessions_anon.RData`: equivalent files for the groups HM2019 and L2019

The objects contained in the COORDINATES files are:

* `{groupyear}_allX`: n_individuals x n_times matrix containing x coordinates (eastings) for each individual at each time point. allX[i,t] gives the x position of individual i at time point t

* `{groupyear}_allY`: same as  above but for y coordinates (northings)

* `{groupyear}_indInfo`: table of information about the individual meerkats (name, DOB, status, domStart, code, dye, sex, shortCode, color). Most relevant columns are code (individual ID code) and status (age/sex/dominance status). The rows in this table correspond to the rows in the allX and allY matrices.

* `{groupyear}_timeLine`: vector of timestamp strings (in UTC) associated with each GPS point

* `{groupyear}_movSummary`: not used

The objects contained in the DATAPRESENCE files are:

* `{groupyear}_audioOn`: matrix of the same dimensions as allX and allY above, indicating whether the audio data for that individual at that time point is available (T if available, F if not). This is important because we have only labeled part of the audio, and also sometimes for the focal-followed meerkats they go out of range of the observer for brief periods, which we exclude.

* `{groupyear}_gpsOn`: same as above, but for GPS data. GPS data may become unavailable if a focal-followed meerkat went out of range of the observer, as described above.


### Scripts

The part contains the following scripts:

* `get_call_response_sequences.R`: code for generating call-response sequences that are later summarized via plots

* `plot_call_response_sequences.R`: uses the saved call-response sequences to produce the plots

* `get_spatiotemporal_call_clustering.R`: code for generating spatiotemporal clustering metric later used in plots

* `plot_spatiotemporal_call_clustering.R`: code for plotting spatiotemporal clustering metric results

### How to run

Once you have installed the required libraries and ensured that you have the data, you are ready to run the analyses.
There are two analyses: an analysis of call and response patterns (as a function of space, time, and individual age/sex/dominance class) and an analysis of the spatial and temporal clustering of calls at the group scale.

#### For the call-response curve analysis:

**Step 1: Generate the call-response sequences that will be used in the plots and save.** This can be generated via the script `get_call_response_sequences.R`. 

* You will need to modify the paths to the data files and directory where you want to save the output on your machine. 

* You may also need to modify two parameters at the top of the script to specify what type of call your are performing the analysis on (by default, the script will run the analysis for close calls). The options are close calls ('cc') or short notes ('sn'). 

* The places where you need to modify are clearly marked at the top of the script with the comment `#YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE`. 

* After this script finishes running, it will save an output file called `callresp_cc_cc_bw0.1.RData` or `callresp_sn_sn_bw0.1.RData` depending on whether you selected close calls (cc) or short notes (sn). This file will be saved in the directory that you specify at the top of the script (`savedir`).

* **NOTE:** This script will take many hours to run on a typical laptop. It is suggested to run it overnight or on a machine you don't mind leaving for awhile! Alternatively, you can use the pre-generated results files provided with the manuscript. We also provide functionality for testing the code by setting the flag `testflag <- T`. This will only generate call-response sequences for a very small subset of the data. It can be used to test the code, however please note that the results will be nonsense! These files will be saved with a `_test` flag at the end of the filename.

**Step 2: Plot call-response curves.** Once you have generated the required output file, you are ready to run the second script: `plot_call_response_sequences.R`. This script will generate the plots from the manuscript. 

* To run this script, you will need to specify the path to your output file (from Part 1) as well as the directory where you have stored the raw data. 

* The places where you need to modify are clearly marked at the top of the script with the comment `#YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE`. 

#### For the spatiotemporal call clustering analysis:

**Step 1: Generate the spatiotemporal call clustering data that will be used in the plots and save.** This can be generated via the script `get_spatiotemporal_call_clustering.R`. 

* You will need to modify the paths to the data files and directory where you want to save the output on your machine. 

* You may also need to modify the parameter `callType` at the top of the script to specify what type of call your are performing the analysis on (by default, the script will run the analysis for close calls). The options are close calls ('cc') or short notes ('sn'). 

* The places where you need to modify are clearly marked at the top of the script with the comment `#YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE`.

* After this script finishes running, it will save output files called `{calltype}_clustering_{groupyear}.RData` for each group-year (HM2017, HM2019, L2019). This file will be saved in the directory that you specify at the top of the script (`savedir`).

* **NOTE:** This script will take several hours to run on a typical laptop. It is suggested to run it on a machine you don't mind leaving for awhile! Alternatively, you can use the pre-generated results files provided with the manuscript. We also provide functionality for testing the code by setting the flag `testflag <- T`. This will only generate clustering metrics for a very small subset of the data. It can be used to test the code, however please note that the results will be nonsense! These files will be saved with a `_test` flag at the end of the filename.

**Step 2: Plot the results**. This can be done using the script `plot_spatiotemporal_call_clustering.R`. 

* To run this script, you will need to specify the path to your output folder (from Part 1), and indciate whether you ran the prior sciprt in test mode (`testflag <- T`) or not (`testflag <- F`). You will also need to specify which call type you would like to plot.

* The places where you need to modify are clearly marked at the top of the script with the comment `#YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE`.

## Call rate vs speed and number of nearby neighbors (Figure 6)

This script produces Figure 6 from precomputed metrics data. You can also compute the metrics from the raw data by setting the flag: 
`compute <- T` at the top.

### Requirements - libraries

See the libraries required at the top of the script.

### Requirements - data

The code requires the following data files to run:

* `metrics_data.RData`

### Scripts

The scripts to produce the first 3 figures in the manuscript are the following:

* `call_rate_vs_speed.R`

### How to run

Run the script to produce the plots. You will need to edit the input and output directories to point to the data file and the output folder where you want the plots to be saved.



