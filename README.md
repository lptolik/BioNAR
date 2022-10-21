
# BioNAR: Biological Network Analysis in R

<!-- badges: start -->
<!-- badges: end -->

We designed BioNAR to support a range of network analysis functionality, 
complementing existing R packages and filling the methodological gaps necessary 
to interrogate biomedical networks with respect to functional and disease 
domains. For that purpose, we do not implement network reconstruction directly 
(unless for synaptic networks), as other tools such as Cytoscape and Network 
Analyst do this already. Rathher, we provide a detailed topologically-based 
network analysis package, enabling the researcher to load networks generated 
from the lab’s own meta-data, thus making the tool as widely applicable and 
flexible as possible. We also provide a synaptic proteome network of our own for 
validation.

## Installation

You can install the development version of BioNAR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lptolik/BioNAR")
```

## Example

Detailed description of the functionality and examples are provided in the
package vignette:
``` r
library(BioNAR)
vignette("BioNAR_overview")
```

