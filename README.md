# Matlab Tool - Design of Brainwaves Detection
#### --- Tool with algorithms to design detection of oscillatory bursts in LFP signals, including burst characterization, synthetic signal generation and analysis of detection

## Overview
#### The work flow has 3 major steps: 1. Characterization, 2. Synthesis, 3. Analysis. Each step does following things:

1.	Select in vivo LFP recording as input. Run the algorithm to fit PSD and obtain SNR and other power spectral properties of the recording signal. Then characterize bursts properties, such as distributions and their fit parameters of amplitude peak (AP), cycle numbers (CN) and burst frequency (BF). Finally save the characterization results.

2.	Load results from previous step. Specify synthesizing parameters. Run an optimization for burst atom amplitude distribution parameters. This step is essential for the synthetic data to reproduce the HT-AP distribution of the real data. There are three types of candidate distribution for the burst atom amplitude: exponential, gamma, lognormal. Exponential is recommended. But users can try out all of them and pick the best fitting. The fitness is indicated by KL divergence or JS divergence (lower divergence value indicates better fitting). Then generate synthetic signal with optimized parameter and see the resulting properties PSD and burst statistics as in step 1. Finally, save the synthesis parameters and data.

3.	Analyze the synthetic data from previous step. Define the ground truth and evaluate detection performance (ROC curve). The relation between detection threshold and true/false positive rates will be shown.


## GUI Instructions

### Step 1: Characterization
#### Project directory:
*	Select the project directory where the results files will be saved in. Use the button to open a dialog box or directly type in the text box.
#### Select Data:
*	"LFP Data" specifies the input data, your in vivo recording. Select variables from your Matlab workspace or type in variable name or expression. The data variable needs to be a vector of LFP signal or a cell array with a vector of LFP segments in each element. The latter format is for cases where some segments with artifacts are removed from the recording.
*	"Sampling Rate" is the recording sampling rate in Hz.
#### Fit PSD:
*	"Oscillation Range" is the frequency range of rhythm you are interested in.
*	"Frequency Range for Fitting" is the range where the algorithm tries to fit a 1/f^a line over the power spectral density. You can adjust the range to find a best fitting for the 1/f^a component. Or you can also use the "Autofit" option to automatically find the range.
*	"Outlier Threshold" is a threshold for determining outliers of the 1/f^a fitting, which is used to find a bump in PSD that indicates significant rhythms. It is measured by median absolute deviation (MAD) of power density in dB/Hz. Try lowering this value if a bump is too small to be found.
*	Click on "Fit" button to run. Fitting results will show in the plot on top right.
#### Characterize Bursts and Fit Probability Distribution:
*	"Number of Histogram Bins" specifies the resolution of the distribution histogram. Adjust this number according to your data size.
*	Click on "Run" button to run. Burst statistics will show in the plot on bottom right. Switch the tabs to view different plots.
#### Save Characterization Data:
*	Type a name for the characterization results and click on "Save" button. A dialog box will open and save the result file under the project directory.
#### Results:
*	A text box where information/message are shown during runtime.

### Step 2: Synthesis
#### Characterization Data:
*	Select the file you saved from step 1 from the dialog box or type in the full path to the file.
*	"Synthetic Sampling Rate" is the sampling rate of the synthetic signal. You can use the same rate as the LFP recording which yields best reproduction of it, or you can specify one. Note that large sampling rate will significantly slow down the process, especially the optimization.
*	"Random Seed" specifies the seed for random number generator.
*	When parameters are set, click on "Load" button before moving to next step.
#### Optimize Amplitude Distribution Parameters:
*	Particle swarm optimization is used. It creates a population of models with different parameters (particle swarm) and they will be attracted toward the minimal cost point in the parameter space after each iteration. The convergence criteria is that the number of stall iterations (consecutive iterations with no lower loss value occurring) reaches 5, or maximal number of iterations (30) is reached.
*	"Synthetic Data Length" is the duration of generated signal. Shorter ones yield shorter runtime while longer ones yield lower variation in resulting distribution meaning more stable solution.
*	"Distribution Types" lists the available candidate distributions for the burst atom amplitude distribution that are used to generate the synthetic bursts. You can select multiple items by holding "ctrl"/"shift"+"mouse left click".
*	"Population" specifies the particle swarm population for each candidate distribution. Select the candidate in the dropdown menu and type in the number in the text box. Clicking on reset button will restore default number. Large population will help convergence but takes longer to run.
*	"Loss Function" is the cost used for the optimization. It can be either KL divergence or JS divergence which indicates the similarity between the HT-AP distributions of the real LFP and synthetic LFP.
*	Click on "Fit" button to run. Note that the optimization can take several minutes to finish. The resulting HT-AP distributions of synthetic data generated by optimized parameters of different candidates will show in the plot on top right, together with the distribution for real data.
#### Generate Synthetic LFP:
*	"Synthetic Data Length" is similar to that in previous step. It can be longer since it runs only once. Longer synthetic data yields lower variance in statistics for latter analysis.
*	"Distribution Type" is the one among the optimized ones you will use to generate the synthetic data.
*	Click on "Run" button to generate. The characterization results of the synthetic data similar to that of the real data in step 1 will show in the plot on bottom right.
#### Save Synthetic Parameters:
*	Save the synthesizing parameters obtained into a file. You can do it right after the optimization step. You can also check the box after you run "Generate Synthetic LFP" step to also save the synthetic LFP data. But the file size will be much larger. Saving the data is not essential for step 3 because the synthetic data can be reproduced solely by the synthesizing parameters. Data generated by only one "Distribution Type" will be saved in one file. You can rerun "Generate Synthetic LFP" with another type and save the data in another file.
*	Type in the case name you want for the file in the text box. The case name will be a suffix following the characterization data file name.
*	Click on "Save" button, a dialog box will show and save the file in the project directory. Note that the synthetic parameter file has to be in the same directory as the characterization data file in order to work properly.
#### Results:
*	A text box where information/message are shown during runtime.
