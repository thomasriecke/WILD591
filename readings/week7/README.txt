Tutorials (Rmarkdown files), R function (R-file) and empirical data (asci-text file) associated with the paper (5 files in total, with 2 supporting jpg files of figures used in tutotials).
Tutorial1.rmd shows estimation bias in simulated dataset (see Box 1 & 2 in paper for details). Format is R-markdown file, which can be opened with for example R Studio.
Tutorial2.rmd illustrates bias due to measurement error and how to account for it (see Box 3 in paper for details). Format is R-markdown file, which can be opened with for example R Studio.
Tutorial3.rmd explains how to analyze the real-world case study of group living benefits in red-winged fairy wrens (see Box 4 in paper for details). Format is R-markdown file, which can be opened with for example R Studio.
Melegans.txt contains the emprical data for red-winged fairy wrens (Malurus elegans) for each of the 108 groups across 9 years. Presented are 698 values for the following variables
- SubjectID: the ID of each group 
- Time:  the year of study with the first year of study 2008 coded as 1
- GroupSize: adult group size during the breeding season
- Survivors: the number of adults in a group that survives till thestart of the next breeding season
- Offspring: the group productivity in terms of number of offspring produced in a year that survives till the next year
- LaggedUnavailable: Boolean denoting whether lagged data is available, thus 1 means  lagged value the lagged value is missing, which is the case for all records in first year of study.
- OffspringLagged:  the value of Offspring in the previous year (Time-1)
- GroupSizeLagged: the value of GroupSize in the previous year (Time-1
There are no further missing values, see description in Box 4 in paper and references therein for details.Format is text file.
simulation_functions.R contains the R functions used in Tutorials 1-3. Format is R code file, which can be opened with for example R Studio or R.
Fig_tutorial1.png a figure used by R markdown in tutorial 1.
Fig_tutorial2.png a figure used by R markdown in tutorial 2.

