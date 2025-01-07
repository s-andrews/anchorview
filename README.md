# anchorview - An R package for visualising Seurat Integration Anchors

This repository contains an R package which is an add-on to the standard [Seurat](https://satijalab.org/seurat/) package for analysing single cell sequencing data.

This package provides tools to construct a visualisation which is useful when performing integration of multiple datasets.  Seurat runs this analysis but doesn't provide a nice way to view the setup of the integration before running it.  This package fills that gap.

![Anchorview Summary](https://raw.githubusercontent.com/s-andrews/anchorview/refs/heads/main/images/anchorview.png)

As well as looking at anchors from individual cells you can also summarise linkages around clusters defined in your two datasets.

![Anchorview Summary](https://raw.githubusercontent.com/s-andrews/anchorview/refs/heads/main/images/anchorview2.png)



## Installation

This package can be installed with the `devtools` or `remotes` R package installers.

### Remotes
```
install.packages("remotes")
remotes::install_github("s-andrews/anchorview")
```

### Devtools
```
install.packages("devtools")
devtools::install_github("s-andrews/anchorview")
```

## Using the package
See the [Anchorview Vignette](https://html-preview.github.io/?url=https://github.com/s-andrews/anchorview/raw/refs/heads/main/inst/doc/anchorview.html) for details of how to use the package.


