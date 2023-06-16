
# RedListEvaluation

<!-- badges: start -->
<!-- badges: end -->

The goal of RedListEvaluation is to use species distribution data to evaluate species against a selection of relevant Red List criteria using the Extent of occurrence (EOO) and Area of occupancy (AOO) metrics. Currently evaluates criteria A2b, A2c, B, D2 based on Red List criteria version 3.1 (see https://www.iucnredlist.org/).

The underlying assumption is that you have distribution data that is representative of your species and area over a time period relevant for a Red List evaluation. 

These methods are based on Brodie et al (2023) "Red List for british seaweeds: evaluating the IUCN methodology for non-standard marine organisms". _Biodiversity and Conservation_.

## Installation

You can install the development version of RedListEvaluation from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cyesson/RedListEvaluation")
```

## Example

This is a basic example which uses some old distribution data to demonstrate the process.

``` r
library(RedListEvaluation)
## basic example code

## load simple dataset of observations
data(Alaria)

## Evaluate all relevant critera for this dataset
Alaria.RLE <- EvaluateAllCriteria(Alaria)

## plot on a map - showing evaluation criteria details (criteria B & D2)
PlotArea(Alaria, Alaria.RLE$B, Alaria.RLE$D2)

## plot trend based on SACFOR abundance observations over the evaluation period
PlotA2b(Alaria.RLE$A2b)

## plot trend in area of occupancy over the evaluation period
PlotA2c(Alaria.RLE$A2c)

```

