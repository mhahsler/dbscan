# DS3: Spatial data with arbitrary shapes

Contains 8000 2-d points, with 6 "natural" looking shapes, all of which
have an sinusoid-like shape that intersects with each cluster. The data
set was originally used as a benchmark data set for the Chameleon
clustering algorithm (Karypis, Han and Kumar, 1999) to illustrate the a
data set containing arbitrarily shaped spatial data surrounded by both
noise and artifacts.

## Format

A data.frame with 8000 observations on the following 2 columns:

- X:

  a numeric vector

- Y:

  a numeric vector

## Source

Obtained from <http://cs.joensuu.fi/sipu/datasets/>

## References

Karypis, George, Eui-Hong Han, and Vipin Kumar (1999). Chameleon:
Hierarchical clustering using dynamic modeling. *Computer* 32(8): 68-75.

## Examples

``` r
data(DS3)
plot(DS3, pch = 20, cex = 0.25)
```
