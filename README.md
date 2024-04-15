# Common neighborhood estimatioj under edge Local Differential Privacy

## Overview

ldp-nb is a C++ project focused on common neighborhood estimation in bipartite graphs with edge local differential privacy.

## Project Structure

The project includes the following files and directories:

- `ldp-nb.cpp`: Implementation of the common neighborhood estimation algorithms. 
- `bigraph.cpp`: Implementation of bipartite graph-related functionality.
- `utility.cpp`: Utility functions used in the project.
- `main.cpp`: Main program entry point.
- `bigraph.h`: Header file for bipartite graph-related functions.
- `utility.h`: Header file for utility functions.
- `ldp-nb.h`: Header file for the common neighborhood estimation algorithms. 
- `makefile`: Build instructions for compiling and linking the project.
- `mt19937ar.h`: Header file for the Mersenne Twister random number generator.

## Build Instructions

To build the project, use the following command:

```bash
make clean && make 
```

## Running the Program

To run the ldp-nb program, use the following command:

```bash
./ldp-nb <epsilon> <data_directory> <num_iterations> <balancing_factor> <sample_size>
```

To compare the performances of all algorithms, set <num_iterations> to 3. Naive and OneR will be run in the first iteration, followed by MultiR-SS and MultiR-DS. 

<balancing_factor> corresponds to kappa, which is used to quantify the degree of imbalance between two vertex degrees in a sampled of vertex pairs.

## Data Format for Bipartite Graphs

This program processes bipartite graph data, which consists of two files: an edge list file (`<datafile>.e`) and a metadata file (`<datafile>.meta`).

### Edge List File (`<datafile>.e`)

The edge list file represents the connections between upper and lower vertices in the bipartite graph. Each line in the file describes an edge between an upper vertex and a lower vertex. The format is as follows:

<upper_vertex> <lower_vertex>

### Metadata File (<datafile>.meta)
The metadata file provides essential information about the bipartite graph:

Upper Vertices Count: The number of upper vertices 
Lower Vertices Count: The number of lower vertices 
Edges Count: The total number of edges 

<upper_vertices_count>
<lower_vertices_count>
<edges_count>

## Usage

1. Run the naive algorithm for 10 rounds with a privacy budget epsilon = 2, on the dataset unicode:  
```bash
./ldp-nb 2 unicode 10 1 btf 
```

2. Run the Multiple-round algorithm for 10 rounds with a privacy budget epsilon = 2, on the dataset unicode:  
```bash
./ldp-nb 2 unicode 10 2 btf 
```