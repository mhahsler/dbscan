# Native sanitizer checks

From the repository root, on Linux with GCC, R, Rcpp, testthat, and the
package's test dependencies installed:

```sh
bash tools/test-native-sanitizers.sh
```

The runner copies the R package into a temporary directory, rebuilds all its
C++ sources with AddressSanitizer and UndefinedBehaviorSanitizer, and runs the
full testthat suite against that installed copy. LeakSanitizer is enabled.
It unloads dbscan before exit to release ANN's package-global KD_TRIVIAL.
Build logs and the installed package stay in the printed temporary directory
for inspection. The working tree's compiled objects are not changed.

The randomized backend tests compare exact kd-tree, linear, and distance-matrix
searches for kNN, fixed-radius neighborhoods, DBSCAN/DBSCAN*, weighted DBSCAN,
LOF, OPTICS, and point density. Explicit query searches are checked against
independent cross-distances. Additional cases cover tied distances, duplicate
observations, zero radius, exact radius boundaries, and errors after native
allocation. The tests allow valid differences in cutoff-tied kNN selections
and ambiguous DBSCAN border assignments.

A standalone C++ harness also injects allocation failures into every stage of
data/query point-storage construction and checks move ownership. It runs with
the same sanitizers and requires no initialized R runtime.
The R driver delivers real SIGINTs during all eight native ANN entry points
(including explicit-query searches) before running the full test suite.

Ordinary tests run 20 random continuous datasets. The sanitizer runner defaults
to 100; override this with `DBSCAN_RANDOMIZED_ITERATIONS`. Seeds are fixed for
reproducibility. All searches use `approx = 0`.
The native backend tests use `skip_on_cran()`. The sanitizer driver explicitly
sets `NOT_CRAN=true` so these tests always run in sanitizer checks. For manual
`R CMD check` runs, set `NOT_CRAN=true` to include them.

The runner uses the C++ compiler configured by R and locates its GCC ASan
runtime. It invokes R's direct executable so ASan is not injected into the R
shell launcher's subprocesses. The C++ runtime is also loaded at startup so
ASan can intercept C++ exception handling. The R executable itself and installed R
dependencies are not rebuilt with sanitizers; package code and inline Rcpp
code compiled into dbscan are instrumented. See R's
[sanitizer guidance](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-the-Address-Sanitizer).

LeakSanitizer needs process-inspection permissions. Containers or sandboxes
that trace or restrict process inspection can prevent its exit-time scan;
run outside that restriction for leak checking. `ASAN_OPTIONS` and
`UBSAN_OPTIONS` can override the defaults when investigating a report.
