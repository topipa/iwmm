# Contributing to IWMM

If you have any ideas for improving this project, or want to contribute something,
feel free to open an issue or pull request.

## Pull request process

*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.
* You can use lintr::lint_package() to check for any styling issues
* You can run styler::style_pkg() to automatically fix styling issues
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat) for testing.
