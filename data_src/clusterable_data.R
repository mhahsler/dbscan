# Create data/clusterable_data.rda from the HDBSCAN notebook file
# clusterable_data.npy:
# https://github.com/scikit-learn-contrib/hdbscan/blob/master/notebooks/clusterable_data.npy
# Convert to CSV first, e.g.:
#   python3 -c "import numpy as np; np.savetxt('clusterable_data.csv',
#     np.load('clusterable_data.npy'), delimiter=',', fmt='%.17g')"

clusterable_data <- as.data.frame(
  read.csv("clusterable_data.csv", header = FALSE, col.names = c("X", "Y")))
save(clusterable_data, file = "../data/clusterable_data.rda", compress = "xz", version = 2)
