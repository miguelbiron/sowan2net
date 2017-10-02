
<!-- README.md is generated from README.Rmd. Please edit that file -->
`sowan2net` : Inference on Network Edge Weights from Sums of Weights At the Nodes

About
-----

The purpose of this package is to find weighted adjacency matrices of a graph or network that match data on the sums of the weights at each node. Since this is an underdetermined problem, the main function of the package is built to explore the solution space.

For a more detailed description of the problem, run `vignette("sowan2net")`.

Installation
------------

``` r
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('miguelbiron/sowan2net')
```

To Do
-----

-   Add working example to readme
-   Improve documentation
