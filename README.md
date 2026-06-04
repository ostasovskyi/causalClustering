# causalClust

`causalClust` provides tools for designing and evaluating cluster-randomized experiments under network interference. The package implements causal clustering procedures for choosing experimental clusters when units interact through a network and the target estimand is a global treatment effect.

The main workflow is:

1. provide a symmetric adjacency matrix for the experimental network;
2. choose a spillover calibration parameter `xi` or a range of calibrations;
3. run a causal clustering algorithm;
4. inspect the bias, variance, and objective components of the selected clustering;
5. optionally compare against baseline clusterings such as epsilon-net, spectral clustering, or Louvain clustering.

## Reference

This package implements methods from:

Viviano, D., Lei, L., Imbens, G., Karrer, B., Schrijvers, O., & Shi, L. (2026). Causal clustering: design of cluster experiments under network interference. Manuscript, May 28, 2026.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("ostasovskyi/causalClustering")
```

Load the package with:

```r
library(causalClust)
```

Most lightweight functionality uses base R and `stats`. Some features require additional packages for graph algorithms, plotting, and sparse matrix operations:

```r
install.packages(c("igraph", "ggplot2", "Matrix"))
```

The package can be installed and used without an SDP solver. The spectral engine provides a lightweight approximation for exploratory computation. It evaluates the same squared-bias objective after rounding, but it does not solve the semidefinite relaxation and does not provide the SDP lower bound or approximation certificate:

```r
fit <- causal_clustering_algorithm1(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 20,
  engine = "spectral",
  methods = "kmeans",
  seed = 123
)
```

The SDP engine implements the paper-style semidefinite relaxation and, when available, returns the SDP-based approximation certificate. It requires the optional package `sdpt3r`.

`sdpt3r` is not currently available from ordinary CRAN installation. Install it from GitHub with:

```r
install.packages("remotes")
remotes::install_github("AdamRahman/sdpt3r")
```

After installing `sdpt3r`, use the SDP engine with:

```r
fit_sdp <- causal_clustering_algorithm1(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 20,
  engine = "sdp",
  methods = "kmeans",
  seed = 123
)
```

## Quick start

The following example uses the spectral approximation engine and does not require an SDP solver.

```r
library(causalClust)

# Simulate a small network with spillovers
sim <- simulate_network_data(
  parameters_graph = list(
    n = 100,
    type_graph = "geometric",
    neighb = 2
  )
)

W <- sim$W

# Run the causal clustering search for one calibration value
fit <- causal_clustering_algorithm1(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 20,
  engine = "spectral",
  objective_type = "squared",
  methods = "kmeans",
  seed = 123
)

fit$selected_k
fit$objective
head(fit$clusters)
fit$components
```

## SDP engine and approximation certificate

When `sdpt3r` is available, use the SDP engine to run the implementation corresponding to the paper's semidefinite relaxation. The returned object includes the SDP lower bound and the approximation certificate when these quantities are available.

```r
fit_sdp <- causal_clustering_algorithm1(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 20,
  engine = "sdp",
  objective_type = "squared",
  methods = "kmeans",
  seed = 123
)

fit_sdp$Gamma_n
fit_sdp$certificate_valid
```

With the SDP engine, the `k_constraint` argument controls how the relaxation is solved:

- `k_constraint = FALSE`: solve one SDP relaxation without optional K-specific constraints and round it for each candidate `K`;
- `k_constraint = TRUE`: solve a separate K-specific SDP relaxation for each candidate `K`, including optional size and box constraints.

## Causal clustering over a calibration range

Use `adaptive_causal_clustering()` or `causal_clustering_algorithm2()` when the spillover calibration is uncertain. These functions use the SDP engine and therefore require `sdpt3r`. The endpoints of `xi_range` must be positive.

```r
fit_adaptive <- adaptive_causal_clustering(
  W = W,
  xi_range = c(0.5, 5),
  min_k = 2,
  max_k = 20,
  engine = "sdp",
  methods = "kmeans",
  seed = 123
)

fit_adaptive$selected_k
fit_adaptive$objective
fit_adaptive$components
```

## Objective components

For a clustering `clusters`, the package computes the objective

```text
xi * variance + bias^2
```

where `variance` is the normalized sum of squared cluster sizes and `bias` is the average fraction of neighbors assigned to different clusters.

```r
W_small <- matrix(
  c(
    0, 1, 1, 0,
    1, 0, 0, 1,
    1, 0, 0, 1,
    0, 1, 1, 0
  ),
  nrow = 4,
  byrow = TRUE
)

clusters <- c(1, 1, 2, 2)

components <- clustering_objective_components(W_small, clusters)
objective <- compute_causal_clustering_objective(W_small, xi = 1, clusters = clusters)

components
objective
all.equal(objective, components$variance + components$bias^2)
```

## Core functions

### Clustering algorithms

| Function | Purpose |
|---|---|
| `causal_clustering_algorithm1()` | Causal clustering algorithm for one calibration value. |
| `causal_clustering_algorithm2()` | Endpoint-regret causal clustering over a calibration range. |
| `adaptive_causal_clustering()` | User-facing wrapper for Algorithm 2. |
| `cluster_epsilon_net()` | Epsilon-net baseline clustering. |
| `cluster_louvain_membership()` | Louvain community detection baseline. |
| `cluster_spectral()` | Spectral clustering baseline. |

### Objective and graph utilities

| Function | Purpose |
|---|---|
| `clustering_objective_components()` | Computes the variance and bias components of a clustering. |
| `compute_causal_clustering_objective()` | Computes `xi * variance + bias^2`. |
| `compute_frontier_objective()` | Wrapper for the Equation (9) frontier objective. |
| `compute_mse_clustering_objective()` | Computes the MSE-style objective using a calibration parameter. |
| `graph_laplacian()` | Computes normalized or unnormalized graph Laplacian. |
| `left_normalized_adjacency()` | Computes the left-normalized adjacency matrix. |
| `available_discretization_methods()` | Lists available rounding/discretization methods. |

### Simulation and evaluation

| Function | Purpose |
|---|---|
| `simulate_network_data()` | Simulates network data with spillovers. |
| `assign_cluster_treatment()` | Assigns treatment at the cluster level. |
| `run_single_network()` | Runs one simulation and compares clustering procedures. |
| `simulate_network_grid()` | Runs a grid of network simulations. |

### Plotting

| Function | Purpose |
|---|---|
| `plot_objective_path()` | Plots objective and selected-cluster paths. |
| `plot_mse_path()` | Plots weighted MSE, bias, variance, and cluster counts. |
| `plot_cluster_count_by_sample_size()` | Plots cluster counts across sample sizes. |
| `plot_objective_grid()` | Plots average objective curves. |
| `plot_mse_grid()` | Plots average weighted-MSE, bias, and variance curves. |

## Input requirements

The main input is an adjacency matrix `W`:

- `W` must be a numeric square matrix;
- `W` should represent an undirected network;
- nonzero entries define the binary adjacency support used by the design criterion;
- diagonal entries are set to zero internally;
- if `W` is not exactly symmetric, it is symmetrized internally as `(W + t(W)) / 2`.

## Calibration parameters

The main tuning parameter is `xi`. In the notation of the causal clustering paper,

```text
xi = (lambda * phibar_n^2 / psibar)^(-1)
```

where:

- `phibar_n` controls the magnitude of spillover effects;
- `psibar` controls the scale of outcome variation;
- `lambda` controls the relative weight placed on squared bias.

Larger `xi` places more weight on the variance component and usually favors more clusters. Smaller `xi` places more weight on reducing cross-cluster exposure and usually favors fewer, larger clusters.

You can also pass `calibration` to `causal_clustering_algorithm1()`, where:

```text
calibration = lambda * phibar_n^2 / psibar
xi = 1 / calibration
```

## Discretization methods

Available methods can be inspected with:

```r
available_discretization_methods()
```

Currently supported methods are:

```r
c(
  "kmeans",
  "hierarchical",
  "spectral_norm_kmeans",
  "spectral_unnorm_kmeans",
  "spectral_norm_hierarchical",
  "spectral_unnorm_hierarchical"
)
```

## Example: comparing clustering methods

The following example uses the spectral engine so that it can run without an SDP solver.

```r
out <- run_single_network(
  seed = 123,
  parameters_graph = list(
    type_graph = "geometric",
    n = 100,
    neighb = 2
  ),
  xi_seq = seq(0.1, 5, by = 0.5),
  min_k = 2,
  max_k = 20,
  include_louvain = TRUE,
  engine = "spectral",
  objective_type = "squared",
  methods = "kmeans"
)

head(out$results)

plots <- plot_objective_path(out$results)
plots$plot_objective
plots$plot_num_clusters
```

## Returned objects

`causal_clustering_algorithm1()` returns a list with entries such as:

- `clusters`: selected cluster labels;
- `selected_k`: selected number of clusters;
- `objective`: selected objective value;
- `components`: bias, variance, number of clusters, and cluster sizes;
- `search_result`: full search output;
- `sdp_lower_bound`: SDP lower bound, when available;
- `Gamma_n`: approximation certificate, when available;
- `certificate_valid`: whether the approximation certificate is available and valid.
