# Example Data for HDBSCAN(e) with Clusters of Various Densities

Contains 118 2-d points with clusters of various densities. Some points
occur more than once. The data is used to show how
[`hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
with `cluster_selection_epsilon` (HDBSCAN(e)) combines the cluster
selection of HDBSCAN\* and DBSCAN\*.

## Format

A data frame with 118 observations on the following 2 variables.

- x:

  a numeric vector

- y:

  a numeric vector

## References

Malzer, C., & Baum, M. (2020). A Hybrid Approach To Hierarchical
Density-based Cluster Selection. In 2020 IEEE International Conference
on Multisensor Fusion and Integration for Intelligent Systems (MFI), pp.
223-228.
[doi:10.1109/MFI49285.2020.9235263](https://doi.org/10.1109/MFI49285.2020.9235263)

## Examples

``` r
data(epsilon_data)
plot(epsilon_data)


## HDBSCAN splits one cluster
hdb <- hdbscan(epsilon_data, minPts = 3)
hullplot(epsilon_data, hdb, main = "HDBSCAN")
#> Warning: Not enough colors. Some colors will be reused.


## HDBSCAN(e) finds all clusters
hdbe <- hdbscan(epsilon_data, minPts = 3, cluster_selection_epsilon = 1)
hullplot(epsilon_data, hdbe, main = "HDBSCAN(e)")
```
