# bmi206-individual-project

This repository contains the data and code required to recreate my analysis for the BMI 206 individual project, based on ["Gene set enrichment analysis for genome-wide DNA methylation data"](https://doi.org/10.1186/s13059-021-02388-x) by Maksimovic et al. (2021). It also contains the specific plots I generated, in the `figures` directory.

### Prerequisites

The `recover-gometh.R` script depends on the following packages: [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [DMRcate](https://bioconductor.org/packages/release/bioc/html/DMRcate.html), [missMethyl](https://bioconductor.org/packages/release/bioc/html/missMethyl.html), [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [tibble](https://cran.r-project.org/web/packages/tibble/index.html), and [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html).

It also expects two large files in the `data` directory from the analysis presented in the paper, which are not in this repository: `GSE110554-data.RData` and `annEPIC.RData`. Both of these are available via [this](https://ucsfonline-my.sharepoint.com/:u:/g/personal/umair_khan_ucsf_edu/EXb4meTgwE5Mp7bgyUKGIkwBWWZCSsBugCLzEBcd7yMdhA?e=xZozaL) link (UCSF login required). Extract the files and place them in the `data` directory before proceeding.

### Running the code

To perform the analysis from start to finish, simply run `recover-gometh.R`. In RStudio, doing so should produce the following three objects:

- `p1` is a plot comparing GOmeth performance through the top 100 ranked terms as the FDR cutoff is reduced from 0.05 to 1 x 10<sup>-4</sup>
- `p2` is a plot charting how the performance of GOmeth changes over the top 100 terms as fewer and fewer CpG sites are sampled from the sets generated with a 0.05 and 1 x 10<sup>-4</sup> FDR cutoff
- `p3` is the same as `p2`, but across only the top 10 terms

### Using pre-computed results

Since `recover-gometh.R` will compute all gene set enrichment analyses from scratch, it will take 10-15 hours to run. To skip this and just regenerate the plots, you can modify the script to use the pre-computed results from my run of `recover-gometh.R`, which are already in the `data` directory.

First, replace line 34 as follows:

```R
dmrGO <- readRDS("data/dmrcate-go.rds")        # old
dmrGO <- readRDS("data/dmrcate-go-new.rds")    # new
```

Then comment out lines 37, 64-84, and 121-154. After making these modifications, the script will generate the same plots in a few seconds or minutes.