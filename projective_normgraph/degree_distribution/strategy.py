# -*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod

import sage.all

from projective_normgraph.util import *


def compute_degrees(ng, degree_distribution, F):
    ng.fe.load_cache()
    for A in F:
        degree_distribution.inc(ng.degree(A))


def compute_restricted_degrees(ng, degree_distribution, F, S):
    ng.fe.load_cache()
    for A in F:
        degree_distribution.inc(ng.restricted_degree(A, S))


class AbstractDegreeDistributionComputationStrategy(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_sample_description(self):
        pass

    @abstractmethod
    def compute(self, norm_graph, pool, r, degree_distribution):
        pass


class RestrictedSubsetDDCS(AbstractDegreeDistributionComputationStrategy):

    def __init__(self, subset, sample_description, subset_nbhd=None):
        self.subset = subset
        self.sample_description = sample_description
        self.subset_nbhd = subset_nbhd

    def get_sample_description(self):
        return self.sample_description

    def compute(self, norm_graph, pool, r, degree_distribution):
        norm_graph.fe.load_cache()
        for (_, NH) in norm_graph.generic_vertex_set_with_neighbourhood_iterator(self.subset, r, self.subset_nbhd):
            degree_distribution.inc(len(NH))


class CanonicalPairDDCS(AbstractDegreeDistributionComputationStrategy):

    def __init__(self, t):
        self.t = t

    def get_sample_description(self):
        return "all_containing_type_{}_canpair".format(self.t)

    def compute(self, norm_graph, pool, r, degree_distribution):
        norm_graph.fe.load_cache()
        pair = norm_graph.canonical_pair(self.t)
        neighbourhood = norm_graph.canonical_pair_neighbourhood(self.t)
        for (_, NH) in norm_graph.generic_vertex_set_with_neighbourhood_iterator(pair, r, neighbourhood):
            degree_distribution.inc(len(NH))


class FullGenericDDCS(AbstractDegreeDistributionComputationStrategy):

    def get_sample_description(self):
        return "all_generic"

    @staticmethod
    def __compute_worker__(pool, norm_graph, r, ddcs, _):
        return norm_graph.degree_distribution(r, ddcs)

    def compute(self, norm_graph, pool, r, degree_distribution):
        if not norm_graph.variant.is_bipartite():
            raise ValueError("This double counting method only works in the bipartite norm graph")
        if r < 2:
            raise ValueError(
                "This double counting method uses canonical pairs and automorphisms to speed up the calculations. "
                "It does not work for r<2.")

        for type_class in norm_graph.type_equivalence_classes():
            norm_graph.canonical_pair_neighbourhood(min(type_class))
        for ((_, _, _, _, type_class_size), result) in pool.run(self.__compute_worker__,
                                                                [(norm_graph, r, CanonicalPairDDCS(min(type_class)),
                                                                  len(type_class)) for type_class in
                                                                 norm_graph.type_equivalence_classes()]
                                                                ):
            result.mul(type_class_size)
            degree_distribution.merge(result)

        degree_distribution.mul(norm_graph.n * len(norm_graph.L.unit_group) / (r * (r - 1)))
        expected_size = sage.all.binomial(norm_graph.L.size(), r) * (len(norm_graph.K.unit_group) ** r)
        assert degree_distribution.size() == expected_size, "%d != %d" % (degree_distribution.size(), expected_size)


class FullDDCS(AbstractDegreeDistributionComputationStrategy):

    def get_sample_description(self):
        return "all"

    def compute(self, norm_graph, pool, r, degree_distribution):
        degree_distribution.merge(norm_graph.degree_distribution(r, FullGenericDDCS()))
        degree_distribution.inc(0, sage.all.binomial(norm_graph.n, r) - degree_distribution.size())


class SimpleRandomSampleDDCS(AbstractDegreeDistributionComputationStrategy):

    def __init__(self, sample_size):
        self.sample_size = sample_size

    def get_sample_description(self):
        return "random_{}".format(self.sample_size)

    def compute(self, norm_graph, pool, r, degree_distribution):
        compute_degrees(norm_graph, degree_distribution,
                        truncate(norm_graph.random_vertex_set_iterator(r), self.sample_size))


class SimpleRandomGenericSampleDDCS(AbstractDegreeDistributionComputationStrategy):

    def __init__(self, sample_size):
        self.sample_size = sample_size

    def get_sample_description(self):
        return "random_generic_{}".format(self.sample_size)

    def compute(self, norm_graph, pool, r, degree_distribution):
        compute_degrees(norm_graph, degree_distribution,
                        truncate(norm_graph.random_generic_vertex_set_iterator(r), self.sample_size))


class RandomGenericSampleDDCS(AbstractDegreeDistributionComputationStrategy):

    def __init__(self, sample_size):
        self.sample_size = sample_size

    def get_sample_description(self):
        return "random_generic_{}".format(self.sample_size)

    @staticmethod
    def __compute_worker__(pool, norm_graph, type, r, count):
        from projective_normgraph.util import truncate
        from projective_normgraph.degree_distribution import DegreeDistribution
        from projective_normgraph.degree_distribution.strategy import compute_restricted_degrees

        degree_distribution = DegreeDistribution(norm_graph, r, None)
        canonical_pair = norm_graph.canonical_pair(type)
        canonical_pair_nbhd = norm_graph.canonical_pair_neighbourhood(type)
        vxset_iterator = norm_graph.random_generic_vertex_set_iterator(r - 2, set(v[0] for v in canonical_pair))
        compute_restricted_degrees(norm_graph, degree_distribution, truncate(vxset_iterator, count),
                                   canonical_pair_nbhd)
        return degree_distribution

    def compute(self, norm_graph, pool, r, degree_distribution):
        if not norm_graph.variant.is_bipartite():
            raise ValueError("This double counting method only works in the bipartite norm graph")
        if r < 2:
            raise ValueError(
                "This double counting method uses canonical pairs and automorphisms to speed up the calculations."
                "It does not work for r<2.")

        # First, randomly sample the type for all r tries and store it in a dictionary C
        type_equivalence_classes = norm_graph.type_equivalence_classes()
        type_class_map = {type: min(type_class) for type_class in type_equivalence_classes for type in type_class}
        type_counts = {min(type_class): 0 for type_class in type_equivalence_classes}
        for _ in xrange(self.sample_size):
            type = type_class_map[norm_graph.K.unit_group.random_element()]
            if type_counts[type] == 0:
                norm_graph.canonical_pair_neighbourhood(type)
            type_counts[type] += 1

        for (_, result) in pool.run(self.__compute_worker__,
                                    [(norm_graph, type, r, type_counts[type]) for type in type_counts]):
            degree_distribution.merge(result)
