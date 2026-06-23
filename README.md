# Genetic diversity in human tissues

Axel Rosendahl Huber, Ferran Muiños, Joan Enric Ramis-Zaldivar, Maria Andrianova, Abel Gonzalez Perez, Núria Lopez-Bigas

## Installation

Clone the repository to a folder in your PC

Open `mutrisk.manuscript.Rrpoj` to start

``` r
# if not installed, install devtools
if (!"pak" %in% rownames(installed.packages())) {
  install.packages("pak")
}

remotes::install_github("im3sanger/dndscv")
remotes::install_github("AxelRosendahlHuber/wintr")
remotes::install_github("AxelRosendahlHuber/mutrisk")

# install all other packages required: 
options(renv.config.pak.enabled = TRUE) # for faster installation
install.packages("renv")
renv::restore()
```

### HEAD
Install MuSical in an environment "python37_musical" using conda or mamba:

Install MuSical in an environment "python37_musical" using conda or mamba: 

``` bash

mamba create -n python37_musical
mamba activate python37_musical
mamba install numpy scipy scikit-learn matplotlib pandas seaborn

# Download musical
cd  /Path/To/MuSiCal
pip install ./MuSiCal
```
