manateessm: R code for fitting the manatee counts state space model of Scolardi et al. (2024)
---------------------------------------------------------------------------------------------

### Contents

Files in this repository:
* manateessm_main.r: the main R script to follow. It loads all required libraries, compiles the Template Model Builder (TMB) C++ scripts, reads in the data .txt file, does all the required pre-processing, and fits the state space model used in Scolardi et al. (2024)
* manateessm_data.txt: simulated data meant to mimic the Sarasota Bay Region data. It has 100 rows, representing time points at a monthly resolution, and includes the following variables:
  - month: time stamp with format YYYY-MM
  - manatee_counts: the count response variable, with NAs for months without any survey
  - number_surveys: the number of aerial surveys run in each month
  - surveyor_experience: a dummy variable representing the experience of the surveyors, with coding 1 = experienced and 0 = less experienced.
  - temperature: daily maximum water temperature.


### Version History

This is manateessm version 0.1. This is the initial release on GitHub.


### References

Scolardi, K. M., Aeberhard, W. H., and Wilkinson, K. A. (2024) Long-term aerial monitoring of Florida manatees, Trichechus manatus latirostris, in a diverse Gulf Coast environment. Submitted


