# my_R_library
Useful generic R functions

## script: calculate_correlations.R
Function: calculate_cors_fast()
- Outputs a flat data.table 
- Calculates correlations, t-statistic, pvalues, standard error
- it is faster in calculating pvalues than what I have found in base R and some available packages such as Hmisc, psych. 
- it can handle one or two different tables

