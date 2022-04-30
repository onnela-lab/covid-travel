# covid-travel
##1. An overview of our folder:

Our GitHub folder here includes all the R code and data for the paper "Framework for assessing and easing global COVID-19 travel restrictions" for the purpose of replicating our paper. This repository has seven sections: Sections 1 to 5 contain the R code used in the paper, Section 6 includes the formatted COVID data sets, and Section 7 contains the CONETTravel R package.

##Some notes on the R codes 

#1. Rcode folders: The R code is included in five separate sections, where Sections 1 to 4 are simulation studies (same order as in the Supplementary Materials of the paper) and Section 5 is the real data analysis.

Most of the R code has been formatted for the High-Performance Cluster of Harvard University (FARS computing at https://www.rc.fas.harvard.edu/). Therefore some modifications may be needed if users want to run the code in their own systems. In the real data analysis section, due to the commercial nature of the flight data from OAG, we cannot make the data publicly available. However, you can plug in your own data to run the model. 

#2. R packages need to be installed: protoABC and CONETTravel

In our ABC estimation procedure, we use the Replenishment Approximate Bayesian Computation (RABC) variant from Drovandi and Pettitt(Biometrics, 2011). The R package for this RABC variant can be installed at https://github.com/AnthonyEbert/protoABC. 

Users also need to install the CONETTravel package at https://github.com/onnela-lab/covid-travel/CONETTravel

## A Final remark:
If you have any difficulty while running the code, please email thle@hsph.harvard.edu or onnela@hsph.harvard.edu.

Thank you for visiting our page. 
