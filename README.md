
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ADARcas

<!-- badges: start -->
<!-- badges: end -->

The goal of `ADARcas` is to provide some functions that allow the user
to measure ADARs activity by using previously developed signatures, or
by computing his own signature if a dataset with its
reconstructed regulatory network is available. ADARs’ activity will be
represented by a Contextual Activity Score (CAS), that is specific for
human neuronal, mouse neuronal or cancer contexts and can be computed
starting from bulk, single-cell RNA-Seq or spatial transcriptomic data.

<p align="center">
<img src="man/figures/image.png" width="700">
</p>

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.19/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.19/biocthis)*.


