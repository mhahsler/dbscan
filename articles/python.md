# Using dbscan from Python

R, the R package **dbscan**, and the Python package **rpy2** need to be
installed. The following example reads the iris data and calls the R
implementation of DBSCAN from Python.

The Python chunks are not run during a normal package build. To execute
them while rendering this vignette, set the environment variable
`DBSCAN_RUN_PYTHON_VIGNETTE=true`. The first run uses **reticulate** to
create a managed Python environment and install **numpy**, **pandas**,
and **rpy2**. The environment is cached for later use, and all Python
chunks run in one shared session.

``` python
import pandas as pd
import numpy as np
import rpy2.robjects as ro

# Prepare data.
iris = pd.read_csv(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data",
    header=None,
    names=["SepalLength", "SepalWidth", "PetalLength", "PetalWidth", "Species"],
)
iris_numeric = iris[["SepalLength", "SepalWidth", "PetalLength", "PetalWidth"]]

# Import the R dbscan package.
from rpy2.robjects import packages
from rpy2.robjects import pandas2ri
dbscan = packages.importr("dbscan")

# Convert the pandas data frame using local conversion rules.
conversion_rules = ro.default_converter + pandas2ri.converter
with conversion_rules.context():
    iris_r = ro.conversion.get_conversion().py2rpy(iris_numeric)

db = dbscan.dbscan(iris_r, eps=0.5, minPts=5)
print(db)
```

Extract the cluster assignment vector as a NumPy array:

``` python
labels = np.asarray(db.rx2("cluster"), dtype=int)
labels
```

Cluster label `0` identifies noise; positive integers identify clusters.
