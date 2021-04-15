# covid-travel
##1. An overview of our folder:

Our GitHub folder here includes all the R codes and data for the paper "Policies for Easing COVID-19 Pandemic Travel Restrictions" with the intended purpose of the research reproducibility. This folder has 7 seperate Sections: Sections 1 to 5 are the R codes used in the paper, Section 6 include the formatted COVID data, and Section 7 is the CONETTravel R packake.

##Some notes on the R codes 

#1. The R codes are included 5 separate sections, with sections 1 to 4 dedicate for the Simulation studies, Section 5 for the Real data analysis.
Most parts of the R codes are in the format of High-Performance Cluster of Harvard University (FARS computing at https://www.rc.fas.harvard.edu/). Therefore some modifications may be needed if users want to run the codes in their own systems. 
In the real data analysis section, due to the commercial reason of the flight data from OAG, we cannot make the data publicly available. However, you can plug in your own data to run the model. 

#2. Some additional notes on some R packages need to be installed:

In our ABC estimation procedure, we use the Replenishment Approximate Bayesian Computation (RABC) variant from Drovandi and Pettitt(Biometrics, 2011). The R package for this RABC variant can be installed at https://github.com/AnthonyEbert/protoABC. 

Users also need to install the CONETTravel package from our folder at https://github.com/onnela-lab/covid-travel/CONETTravel

## A Final remark:
If you have any difficulty while running the code, please email us: thle@hsph.harvard.edu or onnela@hsph.harvard.edu

Thank you for visiting our page. 
