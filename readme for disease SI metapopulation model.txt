Note: you need to download and install Matlab with statistical toolbox to be able to perform the analyses.

One-patch case:

	"onepatchSIdynamics.m" is the matlab code to simulate SI disease model within one patch and save the data output into a csv file.

	"OnpatchSIdynamics.r" is the R code to use the csv file to generate the corresponding figure for the SI disease case within one patch. 

Multiple-patch case:

	"Matlabcode_csvfile_figure1andS1.m": provides the Matlab code to generate the initial csv file used to plot Figure 1 and S1 in this publication. This code is for the 	case when SI model can spread across patches with the migration of host species. Here in the code, I specifically chose 5 patches, but can be modified as N patches. 

	"R Codes in SI model.r" would use the csv file for the SI model spread across 5 patches to plot the figure S1 in this publication. 

