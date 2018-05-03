# Spatial partition-wise/region-wise/segmented/piecewise regression
## Description
* This program automatically detects the change-points (boundaries) of a spatially
varying coefficient linear model, i.e., the coefficients of a linear model varies
among different spatial partitions/regions/segmentations.
* The algorithm used a greedy merging method to iteratively searching for the best
pair of neighboring regions to merge until meeting stopping criteria.
* More details can refer to [post@jay15summer](https://jay15summer.github.io/2018/04/29/spatial-partition-wise-regression.html).
* Note: the algorithm borrows ideas from a paper ['Fast Algorithms for Segmented Regression'](http://proceedings.mlr.press/v48/acharya16.pdf) and a GitHub repository [DataDog/piecewise](https://github.com/DataDog/piecewise).
## Usage
* Main function is: ``[partition_all, partiaion_slt] = spatial_partition_reg(S, X, y, h, v, T)``.
* The inputs of the algorithm are:
   * S: an N-by-2 matrix representing N spatial locations/coordinates;
   * X: an N-by-m matrix of m-dimension independent variables at locations S;
   * Y: an N-by-1 matrix of dependent variable at locations S;
   * h: an integer specifying the initial number of grids at horizontal direction
   (for initial grid-wise segmentation);
   * v: an integer specifying the initial number of grids at vertical direction
   (for initial grid-wise segmentation, the total number of rectangle shape of grid-wise segmentations is h*v);
   * T: a scalar specifying the threshold for stopping the merging.
   If merging a pair of neighboring regions results in the increase-percentage of the fitting error exceeding T, the algorithm will stop,
    i.e., optimal partitioning is found. Usually, 0<T<1, e.g., T = 0.05.
* The outputs of the algorithm are:
   * partition_all: a cell recording the segmentations of every iterations
   (partition_all.segments) and their fitting errors (partition_all.metric);
   * partition_slt: a map representing the selected partitions/segmentations (this program used 'map' to represent partitioning).
   The keys of the map can be obtained via ``cell2mat(keys(partition_slt))``.
   The data [S X y] of every partition can be obtained via ``partition_slt(keys).data``.
   The neighbor of the selected segmentation is ``partition_slt(keys).neighbor``;
   * two plots: 1. percentage-increase of SSE (metric, sum of squared errors); 2. selected final partitioning.
* 'example1.m' gives a simulated example implementing the algorithm.
