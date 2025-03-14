# coFAST

=========================================================================

coFAST is a spatially-aware cell clustering algorithm with cluster significant assessment. It comprises four key modules: spatially-aware cell-gene co-embedding, cell clustering, signature gene identification, and cluster significant assessment.

Check out  our [Package Website](https://feiyoung.github.io/coFAST/index.html) for a more complete description of the methods and analyses. 


Once the coembeddings of  dataset are estimated by coFAST, the package provides functionality for further data exploration, 
analysis, and visualization. Users can:

* Find the signature genes 
* Visuzlize the coembeddings on UMAP space
* Visuzlize the signature genes on UMAP space


# Installation
"coFAST" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}

# For the newest version of coFAST, users can use method 2 for installation.

# Method 2: Install coFAST from Github

if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/coFAST")

# If some dependent packages (such as `scater`) on Bioconductor can not be installed nomrally, use following commands, then run abouve command.
if (!require("BiocManager", quietly = TRUE)) ## install BiocManager
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```



## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 


Tutorials for coFAST method:

* [SRT data analysis for NSCLC](https://feiyoung.github.io/coFAST/articles/CosMx.html)


For the users that don't have set up system properly, the following setup on different systems can be referred.
## Setup on Windows system
First, download [Rtools](https://cran.r-project.org/bin/windows/Rtools/); second, add the Rtools directory to the environment variable.


## Setup on MacOS system
First, install Xcode. Installation about Xcode can be referred [here](https://stackoverflow.com/questions/8291146/xcode-installation-on-mac).


Second, install "gfortran" for compiling C++ and Fortran at [here](https://github.com/fxcoudert/gfortran-for-macOS).


## Setup on Linux  system
If you use conda environment on Linux system and some dependent packages (such as `scater`) can not normally installed, you can search R package at anaconda.org website. We take the `scater` package as example, and its search result is https://anaconda.org/bioconda/bioconductor-scater. Then you can install it in conda environment by following command.
```{Linux}

conda install -c bioconda bioconductor-scater
```
For the user not using conda environment, if  dependent packages (such as `scater`) not normally installed are in Bioconductor, then use the following command to install the dependent packages.
```{Linux}
# install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```
If  dependent packages (such as `DR.SC`) not normally installed are in CRAN, then use the following command to install the dependent packages.
```{Linux}
# install the package on CRAN
install.packages("DR.SC")
```
## Common errors
* When using function `coembedding_umap()`, user may meet the error: "useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE".
Because the `matrixStats` R package remove the argument "useNames=NA" and change the warning to error. Thus, user can install the old version of `matrixStats` by the following code
```{Linux}
# all old versions that are less than 1.1.0  are ok.
# here we take the version 1.1.0 as an example.
remotes::install_version('matrixStats', version='1.1.0') 
```


# Demonstration

For an example of typical coFAST usage, please see our [Package Website](https://feiyoung.github.io/coFAST/index.html) for a demonstration and overview of the functions included in coFAST.

# NEWs
* coFAST version 0.1.0 (2025-03-14)


