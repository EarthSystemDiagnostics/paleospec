# Overview

PaleoSpec is an R package to assist the analysis of variance and power spectra
# License

still private in the ECUS group... later  GPL v3

# Installation 
install the devtools package
Clone the repository to the target directory
```
git clone https://USERNAME@bitbucket.org/ecus/PaleoSpec.git
```
install the package. On the R command line call

```
install("PATH_TO_THE_FIRNR_DIRECTORY/PaleoSpec")
```
now the package is installed and can be loaded with
```
library(PaleoSpec)
```
and the help can be accessed
```
??PaleoSpec
```
# Modification

The source code is in the PATH_TO_THE_FIRNR_DIRECTORY/PaleoSpec/R directory
after modification, the documentation in the files also has to be adapted
(In Emacs-ESS C-c C-o C-o Generate/modify the Roxygen template)

Than compile the documentation again and reinstall the package
```
document("PATH_TO_THE_FIRNR_DIRECTORY/PaleoSpec/")
install("PATH_TO_THE_FIRNR_DIRECTORY/PaleoSpec/")
```

when satisified, commit the changes
(Using ESS and MAGIT: 
'C-C G'
stage the files with 's'
start commit with 'c c'
Write commit message
finish commit with 'C-c C-c'

and finally push them




# Usage
Please refer to the vignette FirnR.pdf for an overview
[Click](https://bitbucket.org/ecus/firnr/raw/b32264a09492d13c0f9d1db7c59ebc1dd54b3b7b/FirnR.pdf)
and the examples in the help file