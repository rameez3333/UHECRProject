# UHECRProject

Data.py: Reads in the Auger and TA data, and calculates the fudge factor for comparing the two.

Ylm.py: Calcualtes the alms and the Cls. The format for the alms is [(0,0),(1,-1),(1,0),(1,1),(2,-2),...]. (0,0) isn't used in chisq calculations. Also, (0,0) is normalized to 1/sqrt(4pi) since the Ylm's are orthonormal.

CR_funcs.py: Contains the exposure function and other miscellaneous CR functions.

Simulated_Events.py: Generates isotropic events with different distributions.

chisq.py: Calculates the chisq between two sets of alms and the p-value. Also can be run to demonstrate the comparison between different data sets.

fig_*: Generates figures.

KMatrix.py contains an alternate method of calculating the alms that is slow.


