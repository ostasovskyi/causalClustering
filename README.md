# causalClust

`causalClust` provides tools for designing and evaluating cluster-randomized experiments under network interference. The package implements causal clustering procedures for choosing experimental clusters when units interact through a network and the target estimand is a global treatment effect.

The main workflow is:

1. provide a symmetric adjacency matrix for the experimental network;
2. choose a spillover calibration `xi`, or provide a calibration grid/range;
3. run the unified causal clustering function;
4. inspect the selected number of clusters, selected rounding method, bias, variance, and objective components;
5. optionally compare against baseline clusterings such as epsilon-net, spectral clustering, or Louvain clustering.

## Reference

This package implements methods from:

Viviano, D., Lei, L., Imbens, G., Karrer, B., Schrijvers, O., & Shi, L. (2026). *Causal clustering: design of cluster experiments under network interference*. Manuscript, May 28, 2026.

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

The recommended engine for the main causal clustering algorithm is the SDP engine. It implements the paper-style semidefinite relaxation and, when available, returns the SDP lower bound and the approximation certificate. The SDP engine requires the optional package `sdpt3r`.

If `install.packages("sdpt3r")` is unavailable in your R setup, install `sdpt3r` from GitHub:

```r
install.packages("remotes")
remotes::install_github("AdamRahman/sdpt3r")
```

On Windows, installing `sdpt3r` from source may require Rtools.

## Quick start with the SDP engine

This example uses the SDP engine. Install `sdpt3r` first, as shown above.

```r
library(causalClust)

# Simulate a small network with spillovers
sim <- simulate_network_data(
  parameters_graph = list(
    n = 60,
    type_graph = "geometric",
    neighb = 2
  )
)

W <- sim$W

# Run causal clustering for one calibration value.
# When several rounding methods are supplied, the package rounds the relaxation
# with each method, evaluates the realized objective, and keeps the best one.
fit <- causal_clustering_algorithm(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 12,
  engine = "sdp",
  methods = available_discretization_methods(),
  seed = 123
)

fit$algorithm
fit$selected_k
fit$selected_method
fit$objective
fit$components
fit$Gamma_n
fit$certificate_valid
head(fit$clusters)
```

The returned `selected_method` reports which discretization method was selected after rounding. The returned `objective` is the realized value of the design criterion for the selected clustering.

## Calibration grids and Algorithm 2

Use the same function when the calibration is uncertain. If `calibration`, `xi`, `calibration_grid`, or `xi_grid` contains more than one value, or if a range is supplied, `causal_clustering_algorithm()` automatically runs the endpoint-regret version of the algorithm.

```r
fit_grid <- causal_clustering_algorithm(
  W = W,
  calibration_grid = c(0.2, 0.5, 1, 2),
  min_k = 2,
  max_k = 12,
  engine = "sdp",
  methods = available_discretization_methods(),
  seed = 123
)

fit_grid$algorithm
fit_grid$xi_range
fit_grid$selected_k
fit_grid$selected_method
fit_grid$objective
fit_grid$rho_values
fit_grid$components
```

Here `calibration = lambda * phibar_n^2 / psibar`, and `xi = 1 / calibration`. For calibration grids, the algorithm uses the endpoint range implied by the grid.

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

### `causal_clustering_algorithm()`

Main user-facing function for causal clustering. It replaces the need to choose manually between the fixed-calibration and calibration-range functions.

| Argument group | Inputs |
|---|---|
| Network | `W`, a numeric square adjacency matrix. Nonzero entries define the binary adjacency support; the diagonal is set to zero internally; nonsymmetric input is symmetrized. |
| Calibration | Provide exactly one of `xi`, `calibration`, `xi_grid`, `xi_range`, `calibration_grid`, or `calibration_range`. A scalar runs the fixed-calibration algorithm. A grid or range runs the endpoint-regret algorithm. |
| Search range | `min_k`, `max_k`, and optionally `n_eig`. These control the candidate number of clusters and the number of eigenvectors used for rounding. |
| Engine | `engine = "sdp"` solves the semidefinite relaxation. `engine = "spectral"` is a lightweight heuristic and is only available for scalar calibrations. |
| Rounding methods | `methods`, a character vector from `available_discretization_methods()`. All requested methods are tried after the relaxation step and the rounded candidate with the smallest realized objective is selected. |
| Optional constraints | `k_constraint`, `gamma_bar`, and `box_constraints` control the optional K-specific size and box constraints for the SDP engine. |
| Reproducibility and storage | `seed`, `try_sign_flip`, `keep_search_matrices`, and `keep_sdp_solutions`. |

Main outputs:

| Output | Meaning |
|---|---|
| `clusters` | Integer cluster labels for the selected design. |
| `selected_k` | Selected number of clusters. |
| `selected_method` | Rounding/discretization method selected after evaluating the realized objective. |
| `objective` | For scalar calibration, the selected value of `xi * variance + bias^2`. For calibration grids/ranges, the selected endpoint-regret criterion. |
| `components` | Bias, variance, number of clusters, and cluster sizes for the selected clustering. |
| `Gamma_n` | SDP approximation certificate when available for the fixed-calibration SDP algorithm. |
| `certificate_valid` | Whether `Gamma_n` is available and valid for the selected call. |
| `sdp_lower_bound` or `endpoint_lower_bounds` | SDP lower bound(s), when available. |
| `search_result`, `solutions`, `rho_values` | Detailed candidate-level output. |

### `clustering_objective_components()`

Computes the two pieces of the design objective for an existing clustering.

Inputs:

- `W`: adjacency matrix;
- `clusters`: cluster label vector of length `nrow(W)`.

Outputs:

- `variance`: `n^{-2} sum_k n_k^2`;
- `bias`: average fraction of neighbors assigned to different clusters;
- `num_clusters`: number of distinct clusters;
- `cluster_sizes`: vector of cluster sizes.

### `compute_causal_clustering_objective()`

Computes the realized scalar objective for an existing clustering.

Inputs:

- `W`: adjacency matrix;
- `xi`: nonnegative calibration value;
- `clusters`: cluster label vector;
- `objective_type`: currently must be `"squared"`.

Output:

- numeric scalar equal to `xi * variance + bias^2`.

## Secondary functions

| Function | Purpose |
|---|---|
| `available_discretization_methods()` | Lists supported rounding/discretization methods. |
| `causal_clustering_algorithm1()` | Backward-compatible fixed-calibration wrapper. |
| `causal_clustering_algorithm2()` | Backward-compatible endpoint-regret wrapper. |
| `adaptive_causal_clustering()` | Backward-compatible wrapper for calibration grids/ranges. |
| `cluster_epsilon_net()` | Epsilon-net baseline clustering. |
| `cluster_louvain_membership()` | Louvain community detection baseline. |
| `cluster_spectral()` | Spectral clustering baseline. |
| `simulate_network_data()` | Simulates network data with spillovers. |
| `assign_cluster_treatment()` | Assigns treatment at the cluster level. |
| `run_single_network()` | Runs one simulation and compares clustering procedures. |
| `simulate_network_grid()` | Runs a grid of network simulations. |
| `plot_objective_path()` | Plots objective and selected-cluster paths. |
| `plot_mse_path()` | Plots weighted MSE, bias, variance, and cluster counts. |

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

You can also pass `calibration`, where:

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

When more than one method is supplied, each method is applied to the relaxation output for each candidate K. The package then computes the realized design objective for every rounded candidate and selects the one with the smallest objective.

## SDP engine and approximation certificate

With `engine = "sdp"`, the returned object includes the SDP lower bound and, when applicable, the approximation certificate:

```r
fit$Gamma_n
fit$certificate_valid
fit$sdp_lower_bound
```

The `k_constraint` argument controls how the relaxation is solved for fixed-calibration calls:

- `k_constraint = FALSE`: solve one SDP relaxation without optional K-specific constraints and round it for each candidate `K`;
- `k_constraint = TRUE`: solve a separate K-specific SDP relaxation for each candidate `K`, including optional size and box constraints.

## Lightweight spectral approximation

The package can also be used without an SDP solver. The spectral engine provides a lightweight approximation for exploratory computation. It evaluates the same squared-bias objective after rounding, but it does not solve the semidefinite relaxation and does not provide the SDP lower bound or approximation certificate.

```r
fit_spectral <- causal_clustering_algorithm(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 12,
  engine = "spectral",
  methods = available_discretization_methods(),
  seed = 123
)

fit_spectral$selected_k
fit_spectral$selected_method
fit_spectral$objective
fit_spectral$components
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
  methods = available_discretization_methods()
)

head(out$results)

plots <- plot_objective_path(out$results)
plots$plot_objective
plots$plot_num_clusters
```