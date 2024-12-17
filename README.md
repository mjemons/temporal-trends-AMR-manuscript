# Temporal trends in antibiotic resistance in Europe, 1998-2019

This repository organises the code for the [preprint](https://medrxiv.org/cgi/content/short/2023.09.27.23296241v1)

## Folder contents

- `Masterscript.R`: overall script which calls all the individual scripts. The idea here is to make the workflow clear to an external person using the scripts.

- `/code`: folder containing code to generate temporal fits and correlations.
code/preprocessing: code that generates the aggregated datasets from the raw data (not shared).
code/obsolete: older analyses no longer relevant to the paper, but kept for reference (for now)

- `/analysis`: folder containing code to further analyse data sets resulting from model fitting and correlation calculations (i.e. generate figures, tables etc).

- `supporting_analysis`: folder containing supporting analyses - either as part of our thinking or to go in supplementary materials

- `package management`: Packages are managed using `renv`. All package versions are specified in the `renv.lock` and can easily be installed with this framework in your system. `renv` should automatically install itself and the required packages. In case this does not work, do the following: `install.packages("renv")` and then `renv::restore()`.
Futhermore, we provide the `R sessionInfo()` that was used to generate the results. 

- `author contribution`: This package shows the state of the repository at publication. The history was truncated and therefore the contributions are not reflective of the project. The contributors to the code are (in alphabetical order) [François Blanquart](https://sites.google.com/site/francoisblanquart/), [Martin Emons](https://www.mls.uzh.ch/en/research/robinson/groupmembers/martin-emons.html) and [Sonja Lehtinen](https://sites.google.com/view/sonjalehtinen).