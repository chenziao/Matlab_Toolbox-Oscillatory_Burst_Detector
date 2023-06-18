# Matlab Tool - Design of Brainwaves Detection
#### --- Tool with algorithms to design detection of oscillatory bursts in LFP signals, including burst characterization, synthetic signal generation and analysis of detection

## Overview
#### The work flow has 3 major steps: 1. Characterization, 2. Synthesis, 3. Analysis, as follows:

1.	Select local field potential (LFP) recording as input. Run the algorithm to fit power spectral density (PSD), decompose the PSD into signal and background components. Then characterize bursts properties, including statistics of amplitude peaks, number of cycles and burst frequency. Finally, save the characterization results.

2.	Load results from the previous step. Specify parameters for synthesizing the background component and oscillatory signal component. Then perform an optimization for amplitude peak distribution parameters of burst atoms to be generated in the signal component. This step is essential for the synthetic data to reproduce the amplitude distribution of the observed data. There are three candidate types of parametric distribution for the burst atom amplitude: gamma, lognormal and exponential. Gamma distribution works well in general. The goodness of fit is indicated by Kullback–Leibler divergence (D<sub>KL</sub>) or Jensen–Shannon divergence (D<sub>JS</sub>). Lower divergence value indicates better fit. Then synthetic LFP is generated with optimized parameters and the resulting properties including PSD and burst statistics are visualized as in step 1. Finally, save the synthesis parameters and data.

3.	Analyze the synthetic data generated from the previous step for evaluating detection. Define the ground truth with lower and upper bounds in the synthetic signal trace and evaluate detection performance on the composite trace using a receiver operating characteristic (ROC) curve. The relation between detection threshold and true/false positive rates will be shown.


## Flow chart for step 1 and 2

<img src="https://raw.githubusercontent.com/chenziao/Matlab_Tool-Design_of_Brainwaves_Detection/main/image/FlowChart.png" width="625" height="498">


## User Mannual

See the document [User Manual.pdf](User%20Manual.pdf).
