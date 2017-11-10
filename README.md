# SageMath package for projective norm graphs

This SageMath package provides utilities to do computations with projective norm graphs as introduced by [1] and [2]. See there for definitions.
    
## Installation
1. Ensure you have [SageMath](https://www.sagemath.org/) installed. This package was tested on SageMath 8.1. Newer versions probably work as well. Older versions may or may not work.

2. Download the source from this repository.

    $ git clone https://github.com/TomasBayer/projective-normgraph.git

3. Change to the root directory and run

    $ sage -pip install --no-index -v .

## Usage
When installed as described above, this package can be used by the SageMath instance it was installed with like this:

    sage: import projective_normgraph

Or, using any other python import statements for that matter.

    sage: from projective_normgraph import NormGraph as NG

This will work from everywhere. You don't need to work in the root directory of this repository.

## Example
    import projective_normgraph.config as config

    config.set_cache_dir("~/normgraph_cache/")
    config.set_results_dir("~/normgraph_results/")

    # Some operation allow calculations to be carried out simultaneously
    # by multiple processing units. This sets the maximal number of 
    # processes that may be started simultaneously.
    config.set_max_processes(2)

    # Choose a degree distribution computation strategy (DDCS)
    from projective_normgraph.degree_distribution.strategy import *

    # This DDCS computes the distribution over all vertex sets taking 
    # advantage of symmetries to speed up calculation.
    ddcs = FullDDCS() 

    # This DDCS chooses a sample of 500 random generic vertex sets 
    #ddcs = RandomGenericSampleDDCS(500)

    from projective_normgraph import BipartiteNormGraph as BNG
    # iterate over all bipartite projective norm graphs of dimension 4.
    for A in BNG.it(4):
        print A

        # Compute a distribution of degrees of 4-vertex sets
        # using the previously defined strategy
        degree_distribution = A.degree_distribution(4, ddcs)
        print degree_distribution.print_statistics()

## References

- [1]  Kollár, J., Rónyai, L. and Szabó, T. "Norm-graphs and Bipartite Turán Numbers", Combinatorica 16 (1996), no. 3, 399-406
- [2]  Alon, N., Rónyai, L., and Szabó, T. "Norm-graphs: variations and applications", J. Combin. Theory Ser. B 76 (1999), pp. 280–290.
