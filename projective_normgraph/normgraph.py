# -*- coding: utf-8 -*-
import logging
import os
from enum import Enum
from itertools import combinations, ifilter, imap, dropwhile

import sage.all

from projective_normgraph import config
from projective_normgraph.degree_distribution import DegreeDistribution
from projective_normgraph.finitefield import FiniteFieldExtension
from projective_normgraph.processing import ProcessingPool
from projective_normgraph.random import RandomGraph, RandomBipartiteGraph
from projective_normgraph.util import all_distinct, chunk, prime_powers, produce, unique_element, skip_until

logger = logging.getLogger(__name__)


class NormGraphVariant(Enum):
    NG = (0, '+', 'NormGraph')
    MNG = (1, '-', 'NormGraph-')
    BNG = (2, 'b+', 'BipartiteNormGraph+')
    BMNG = (3, 'b-', 'BipartiteNormGraph')

    def __str__(self):
        return self.value[2]

    def is_bipartite(self):
        return self.value[0] > 1

    def is_minus(self):
        return self.value[0] % 2 == 1

    @staticmethod
    def of(variant):
        if isinstance(variant, NormGraphVariant):
            return variant
        for v in NormGraphVariant:
            if variant in v.value:
                return v
        raise ValueError('Illegal variant name: {}'.format(variant))


class NormGraphCanonicalPairNeighbourhoodCache(object):

    # TODO: Do not double store BipartiteNormGraph+ and BipartiteNormGraph-

    def __init__(self, norm_graph):
        self.norm_graph = norm_graph
        self._file_format = "CanonicalNeighbourhood_{!r}_{}.sobj"
        self._cache = {}

    def __get_nbhd_cache_file(self, type):
        cache_dir = config.get_cache_dir()
        if not os.path.isdir(cache_dir):
            raise ValueError("Cache directory does not exist.")
        return os.path.join(cache_dir, self._file_format.format(self.norm_graph, type))

    @staticmethod
    def __build_cache_worker(_, norm_graph, A, v, Y):
        partition = [[] for _ in norm_graph.K.unit_group]
        norm_graph.fe.load_cache()
        for X in A:
            x = norm_graph.edgecomplete_right_vertex(v, X)
            if x is not None:
                w = (X, x)
                y = norm_graph.edgecomplete_left_vertex(Y, w)
                if y is not None:
                    partition[y].append(w)
        return partition

    def __build_cache(self):
        logger.debug("Building canonical neighbourhoods cache...")

        self.norm_graph.fe.load_cache()

        canonical_pair_types = [min(T) for T in self.norm_graph.type_equivalence_classes()]
        v = self.norm_graph.canonical_vertex()
        Y = self.norm_graph.L.one()

        self._cache = {t: [] for t in canonical_pair_types}
        pool = ProcessingPool(max_processes=config.get_max_processes())
        for (_, ret) in pool.run(NormGraphCanonicalPairNeighbourhoodCache.__build_cache_worker,
                                 [(self, A, v, Y) for A in
                                  chunk(self.norm_graph.L,
                                        sage.all.ceil(self.norm_graph.L.size() / pool.max_processes))]):
            for canonical_pair_type in canonical_pair_types:
                self._cache[canonical_pair_type] += ret[canonical_pair_type]

        logger.debug("Saving cache...")
        for canonical_pair_type in canonical_pair_types:
            sage.all.save(self._cache[canonical_pair_type], self.__get_nbhd_cache_file(canonical_pair_type))
        logger.debug("Cache saved.")

    def __load_nbhd_cache(self, canonical_pair_type):
        if canonical_pair_type not in self._cache:
            cache_file = self.__get_nbhd_cache_file(canonical_pair_type)
            if os.path.isfile(cache_file):
                logger.debug("Loading cached canonical neighbourhood of type {} in {!r}".
                             format(canonical_pair_type, self))
                self._cache[canonical_pair_type] = sage.all.load(cache_file)
                logger.debug("Cache loaded")
            else:
                self.__build_cache()

    def get(self, canonical_pair_type=1):
        """Returns the neighbourhood of the canonical pair of type t."""
        if not self.norm_graph.variant.is_bipartite():
            raise ValueError("Canonical pairs only exist in the bipartite norm graph!")
        if canonical_pair_type in self._cache:
            return self._cache[canonical_pair_type]
        else:
            canonical_pair_type = min(self.norm_graph.type_equivalence_class(canonical_pair_type))
            self.__load_nbhd_cache(canonical_pair_type)
            return self._cache[canonical_pair_type]


class NormGraph(object):
    """This class represents a projective norm graph or one of its variants.

    Vertices of the norm graph are represented as pairs (X,x) where X is an element from the finite field with q^(t-1)
    elements and x is an element from the unit group of the finite field with q elements.
    
    Vertices of the bipartite variant are also represented as such pairs, a vertex itself does not tell in which part
    it is contained in. This always depends on the context it occurs in. When vertices or vertex sets occur in function
    arguments or return values, it will usually either be clear from the context or irrelevant in which part they are
    contained in. Whenever this is not the case, it will be stated explicitly how function arguments or return values
    should be interpreted.
    
    To that end, the two parts of the bipartition are called "left part" and "right part" and accordingly, each vertex
    is a "left vertex" or a "right vertex". Wherever vertex sets occur, they are implicitly considered to be contained
    entirely in one of the two parts. Accordingly, a vertex set is either a "left vertex set" or a "right vertex set".
    """

    def __init__(self, q, t, variant=NormGraphVariant.NG):
        self.q = q
        self.t = t
        self.variant = NormGraphVariant.of(variant)

        if t < 2:
            raise ValueError("t should be at least 2.")

        F = sage.all.Integer(q).factor()
        if len(F) != 1:
            raise ValueError("q={} is not a prime power.".format(q))

        self.p = F[0][0]
        self.k = F[0][1]
        self.d = self.k * (self.t - 1)

        # NormGraph- only exists ifb subtraction is commutative.
        if self.variant == NormGraphVariant.MNG and t % 2 == 0 and self.p > 2:
            raise ValueError("The norm graph NormGraph-(q,t) is only defined for odd t or q which are a power of 2.")

        self.fe = FiniteFieldExtension(self.p, self.d, self.k)
        self.L = self.fe.L
        self.K = self.fe.K

        self.op_norm = self.fe.difference_norm if self.variant.is_minus() else self.fe.sum_norm

        # The number of vertices in one part
        self.n = self.L.size() * (self.q - 1)

        self._canonical_pair_neighbourhood_cache = NormGraphCanonicalPairNeighbourhoodCache(self)

    @staticmethod
    def it(t, variant=0, min_q=3, max_q=sage.all.infinity, primes_only=False):
        """Returns a generator that yields all t-dimensional norm graphs of the specified variant with q >= min and q
        <= max in ascending order. If primes_only is True, non-prime q are omitted. """
        for q in sage.all.primes(min_q, max_q) if primes_only else prime_powers(min_q, max_q):
            yield NormGraph(q, t, variant)

    def __getstate__(self):
        return self.q, self.t, self.variant

    def __setstate__(self, state):
        self.__init__(*state)

    def __repr__(self):
        """Returns the symbolic name of the norm graph."""
        return "{}({},{})".format(self.variant, self.q, self.t)

    def clone(self):
        return NormGraph(self.q, self.t, self.variant)

    # Graph parameters.

    def num_vertices(self):
        """Returns the number of vertices."""
        return self.n * (1 + int(self.variant.is_bipartite()))

    def num_edges(self):
        """Returns the number of edges."""
        return (self.n * (self.L.size() - 1) - self.num_loop_vertices()) // (1 + int(not self.variant.is_bipartite()))

    def has_loop_vertices(self):
        return self.variant == NormGraphVariant.NG and self.p != 2

    def num_loop_vertices(self):
        """Returns the number of loop vertices."""
        return (self.L.size() - 1) if self.has_loop_vertices() else 0

    def average_degree(self):
        """Returns the average degree."""
        return (self.L.size() - 1) - self.num_loop_vertices() / self.n

    def density(self):
        """Returns the density."""
        return self.num_edges() / sage.all.binomial(self.num_vertices(), 2)

    def num_vertex_sets(self, r):
        """For non-bipartite norm graphs, this returns the number of r-vertex sets. For bipartite norm graphs,
        it returns the number of r-vertex sets contained in one part. """
        return sage.all.binomial(self.n, r)

    def num_generic_vertex_sets(self, r):
        """For non-bipartite norm graphs, this returns the number of generic r-vertex sets. For bipartite norm
        graphs, it returns the number of generic r-vertex sets contained in one part. """
        return sage.all.binomial(self.L.size(), r) * (self.q - 1) ** r

    # TODO: diameter, clique number, independence number, distance

    def max_degree(self, k):
        """Returns the maximal degree of k-vertex sets non-computationally, if possible"""
        if k < 3:
            return (self.L.size() - 1) // (self.K.size() ** k)
        if k == 3:
            if self.t == 4:
                return 2 * self.q + 2 - int((self.q + 1) % 3)
        raise ValueError("The maximal degree of {}-vertex sets is unknown for this graph.".format(k))

    # Vertex and vertex set properties

    def is_generic(self, A):
        """Checks whether A is a generic vertex set."""
        return all_distinct(imap(lambda v: v[0], A))

    def is_aligned(self, A):
        """Checks whether A is an aligned vertex set."""
        return unique_element(imap(lambda v: v[1], A)) is not None

    def is_loop_vertex(self, v):
        """Checks whether the vertex v is a loop vertex."""
        if not self.has_loop_vertices():
            return False
        return not self.L.is_zero(v[0]) and v[0] + self.fe.zech_logarithm_norm(0) == 2 * v[1]

    def pair_type(self, v, w):
        """Returns the type of the pair (v,w)."""
        return w[1] - v[1]

    def triple_type(self, v, w, z):
        """Returns the type of the generic triple (v,w,z)."""
        c1 = v[1] - z[1] + self.fe.difference_norm(w[0], z[0]) - self.fe.difference_norm(v[0], w[0])
        c2 = w[1] - z[1] + self.fe.difference_norm(v[0], z[0]) - self.fe.difference_norm(v[0], w[0])
        return c1, c2

    # Canonical vertex, vertex-sets, neighbourhood

    def canonical_vertex(self):
        """Returns the canonical vertex."""
        return self.L.zero(), self.K.one()

    def canonical_pair(self, t=1):
        """Returns the canonical pair of type t."""
        return self.canonical_vertex(), (self.L.one(), self.K.product(self.K.one(), t))

    def canonical_triple(self, c1=0, c2=0):
        """Returns the canonical triple of type (c1, c2)."""
        if self.q < 3:
            raise ValueError()
        A3 = self.L.generator()
        a3 = self.K.quotient(self.fe.difference_norm(A3, self.L.one()), c1)
        t = self.K.quotient(self.K.product(c2, a3), self.fe.norm(A3))
        return self.canonical_pair(t) + ((A3, a3),)

    def canonical_pair_neighbourhood(self, t=1):
        """Returns the neighbourhood of the canonical pair of type t."""
        return self._canonical_pair_neighbourhood_cache.get(t)

    # Random vertex and vertex set sampling

    def random_vertex(self):
        """Picks a vertex uniformly at random and returns it."""
        return self.L.random_element(), self.K.unit_group.random_element()

    def random_vertex_set(self, k, E=None):
        """Among all k-vertex-sets which are disjoint from E, picks one uniformly at random and returns it."""
        if E is None:
            E = set()
        S = [None for _ in xrange(k)]
        E = list(E)
        for i in xrange(k):
            S[i] = ifilter(lambda x: x not in E, produce(self.random_vertex)).next()
            E.append(S[i][0])
        return tuple(S)

    def random_generic_vertex_set(self, k, E=None):
        """Among all generic k-vertex-sets whose set of first coordinates is disjoint from E, picks one uniformly at
        random and returns it. """
        if E is None:
            E = set()
        S = [None for _ in xrange(k)]
        E = list(E)
        for i in xrange(k):
            S[i] = (ifilter(lambda x: x not in E, produce(self.L.random_element)).next(),
                    self.K.unit_group.random_element())
            E.append(S[i][0])
        return tuple(S)

    # Vertex and vertex sets iterators

    def vertex_iterator(self):
        """Yields all vertices of one partition (in the bipartite case) or all vertices (in the non-bipartite case)."""
        for X in self.fe.L:
            for x in self.K.unit_group:
                yield (X, x)

    def vertex_set_iterator(self, k, skip_until_element=None):
        """Yields all k-element vertex sets within one partition (in the bipartite case) or all k-element vertex sets
        (in the non-bipartite case). """
        if k == 0:
            yield ()
        else:
            for v in self.vertex_iterator() if skip_until_element is None else skip_until(self.vertex_iterator(),
                                                                                          skip_until_element):
                for S in self.vertex_set_iterator(k - 1, v):
                    yield S + (v,)

    def generic_vertex_set_iterator(self, k, skip_until_element=None):
        """Yields all generic k-element vertex sets within one partition (in the bipartite case) or all generic
        k-element vertex sets (in the non-bipartite case). """
        if k == 0:
            yield ()
        else:
            for X in self.L if skip_until_element is None else skip_until(self.L, skip_until_element):
                for x in self.K.unit_group:
                    for S in self.generic_vertex_set_iterator(k - 1, X):
                        yield S + ((X, x),)

    # Random vertex set iterators  

    def random_vertex_set_iterator(self, k, E=None):
        """Endlessly draws samples from all k-vertex-sets which are disjoint from E uniformly at random and yields
        them. """
        if E is None:
            E = set()
        if len(E) + k > self.n:
            raise ValueError()
        while True:
            yield self.random_vertex_set(k, E)

    def random_generic_vertex_set_iterator(self, k, E=None):
        """Endlessly draws samples from all generic k-vertex-sets whose set of first coordinates is disjoint from E
        uniformly at random and yields them. """
        if E is None:
            E = set()
        if k > self.L.size() - len(E):
            raise ValueError("Cannot produce generic {}-vertex sets "
                "if there are only {}-{}={} first coordinates to choose from".format(
                k, len(self.L), len(E), len(self.L) - len(E)))
        while True:
            yield self.random_generic_vertex_set(k, E)

    # Neighbourhood iterator

    def generic_vertex_set_with_neighbourhood_iterator(self, S, k, N=None, skip_until_element=None):
        """For each generic k-element vertex set containing S, yields a pair containing the vertex set and its
        neighbourhood. If the neighbourhood of S is already known, it may be provided as N to reduce computation
        time. """
        if k < len(S):
            return
        if N is None:
            N = self.neighbourhood(S)
        if k == len(S):
            yield (S, N)
        else:
            forbidden = map(lambda v: v[0], S)
            for X in self.L if skip_until_element is None else dropwhile(lambda x: x != skip_until_element, self.L):
                if X not in forbidden:
                    for (v, NN) in self.neighbourhood_partition(X, N):
                        for C in self.generic_vertex_set_with_neighbourhood_iterator(S + (v,), k, NN, X):
                            yield C

    # Degree distribution

    def degree_distribution(self, r, ddcs, cache_load=True, cache_save=True):
        """Returns the distribution of the degrees of r-element vertex sets according to the provided
        DegreeDistributionComputationStrategy ddcs. If the result of this method is available in the cache and
        cache_load is set to True, this returns the result from the cache, otherwise it will be computed. If the
        result is computed, it will be stored in the cache prior to being returned unless cache_save is set to False.
        """
        results_dir = config.get_results_dir()
        if cache_load:
            try:
                return next(DegreeDistribution.load(self, r, ddcs.get_sample_description(), results_dir))
            except StopIteration:
                pass
        degree_distribution = DegreeDistribution(self, r, ddcs.get_sample_description())
        pool = ProcessingPool(max_processes=config.get_max_processes())
        ddcs.compute(self, pool, r, degree_distribution)
        if cache_save:
            degree_distribution.save(ddcs.get_sample_description(), results_dir)
        return degree_distribution

    # Neighbourhood computation

    def neighbourhood(self, A):
        """Returns the set of common neighbours of the vertex set A. (If this norm graph is bipartite,
        A is considered to be a left vertex set.) """
        if not self.is_generic(A):
            return []
        if len(A) == 0:
            return list(self.vertex_iterator())
        neighbours = []
        for X in self.L:
            x = self.starcomplete_right_vertex(A, X)
            if x is not None:
                neighbours.append((X, x))
        return neighbours

    def degree(self, A):
        """Returns the common degree of the vertex set A. (If this norm graph is bipartite, A is considered to be a
        left vertex set.) """
        if not self.is_generic(A):
            return 0
        else:
            return sum(1 for X in self.L if self.starcomplete_right_vertex(A, X) is not None)

    def degree_analytically(self, A):
        """Returns the common degree of the vertex set A non-computationally if possible. (If this norm graph is
        bipartite, A is considered to be a left vertex set.) """
        if len(A) == 1:
            return self.L.size() - 1 - self.is_loop_vertex(A[0])
        if len(A) == 2:
            if A[0] == A[1]:
                raise ValueError("The vertices in A should be distinct, but {} appears twice".format(A[0]))
            if not self.is_generic(A):
                return 0
            vertex_count = (self.L.size() - 1) // (self.K.size() - 1) - int(self.pair_type(A[0], A[1]) == self.K.one())
            if self.is_edge(A[0], A[1]):
                for v in A:
                    vertex_count -= int(self.is_loop_vertex(v))
            return vertex_count
        if len(A) > 2:
            raise ValueError("The degree of {}-vertex sets is unknown (yet)".format(len(A)))
        # TODO: It is not.

    def restricted_degree(self, A, S):
        """Returns the number of common neighbours of the vertex set A within the vertex set S. (If this norm graph
        is bipartite, A is considered to be a left vertex set and S is considered to be a right vertex set.) """
        if not self.is_generic(A):
            return 0
        else:
            return sum(1 for w in S if self.is_rightstar(A, w))

    def neighbourhood_partition(self, X, S):
        """Yields pairs (v, N) where v is a left vertex with first coordinate X and N is the neighbourhood of v
        intersected with S """
        A = {x: [] for x in self.K.unit_group}
        for v in S:
            x = self.edgecomplete_left_vertex(X, v)
            if x is not None:
                A[x].append(v)
        for x in A:
            yield ((X, x), A[x])

    # other

    def sage_graph(self):
        """Returns a sage graph object which is isomorphic to this norm graph."""
        self.fe.load_cache()
        edges = []
        is_bipartite = self.variant.is_bipartite()
        for X in self.L:
            for Y in self.L:
                a = self.op_norm(X, Y)
                if not self.K.is_zero(a):
                    for c in self.K.unit_group:
                        d = a - c
                        if X != Y or c < d or is_bipartite:
                            edges.append(((X, c, False), (Y, d, is_bipartite)))
                if X == Y and not is_bipartite:
                    break
        return sage.all.Graph(edges)

    def automorphism_group_order(self):
        """Returns the order of the automorphism group"""
        if self.q == 2:
            return sage.all.factorial(self.L.size()) * (1 + int(self.variant.is_bipartite()))
        elif self.variant.is_bipartite():
            return (self.L.size() - 1) * self.d * self.L.size() * (self.K.size() - 1) * 2
        elif self.variant.is_minus() or self.p == 2:
            return (self.L.size() - 1) * self.d * self.L.size()
        else:
            return (self.L.size() - 1) * self.d

    def automorphism_group(self):
        pass

    def associated_random_graph(self):
        if self.variant.is_bipartite():
            return RandomGraph(self.n, self.density())
        else:
            return RandomBipartiteGraph(self.n, self.n, self.num_edges() / (self.num_vertices() ** 2))

    # adjacency checks

    def is_edge(self, v, w):
        """Decides whether the vertices v and w are adjacent. (If this norm graph is bipartite, v is considered to be
        a left vertex and w is considered to be a right vertex.) """
        return self.op_norm(v[0], w[0]) == (v[1] + w[1]) and (self.variant.is_bipartite() or v != w)

    def is_clique(self, A):
        """Decides whether the vertex set A forms a clique. (This method is only available if this norm graph is not
        bipartite.) """
        if self.variant.is_bipartite():
            raise ValueError()
        return all(self.is_edge(v, w) for (v, w) in combinations(A, 2))

    def is_independent(self, A):
        """Decides whether the vertex set A is independent. (This method is only available if this norm graph is not
        bipartite.) """
        if self.variant.is_bipartite():
            raise ValueError()
        return not any(self.is_edge(v, w) for (v, w) in combinations(A, 2))

    def is_complete(self, A, B):
        """Decides whether the vertex set A and B forms a complete bipartite subgraph. (If this norm graph is
        bipartite, A is considered to be a left vertex set and B is considered to be a right vertex set.) """
        return all(self.is_edge(v, w) for v in A for w in B)

    def is_leftstar(self, v, B):
        """Decides whether the vertex v is adjacent to all vertices in the vertex set B. (If this norm graph is
        bipartite, v is considered to be a left vertex and B is considered to be a right vertex set. Note,
        that if this norm graph is not bipartite, then leftstar(v, B) == is_rightstar(B, v).) """
        return all(self.is_edge(v, w) for w in B)

    def is_rightstar(self, A, w):
        """Decides whether the vertex w is adjacent to all vertices in the vertex set A. (If this norm graph is
        bipartite, A is considered to be a left vertex set and w is considered to be a right vertex. Note,
        that if this norm graph is not bipartite, then rightstar(A, w) == is_leftstar(w, A).) """
        return all(self.is_edge(v, w) for v in A)

    def edgecomplete_right_vertex(self, v, Y):
        """Given a vertex v and Y in L, returns the y in K, such that v is adjacent to (Y,y), if such y exists (In
        this case it is unique.). Otherwise returns None. (If this norm graph is bipartite, v is considered to be a
        left vertex and (Y,y) is considered to be a right vertex. Note, that if this norm graph is not bipartite,
        then edgecomplete_right_vertex(v, Y) == edgecomplete_left_vertex(Y, v).) """
        h = self.op_norm(v[0], Y)
        if self.K.is_zero(h):
            return None
        h -= v[1]
        if self.variant == NormGraphVariant.NG and v == (Y, h):
            return None
        return h

    def edgecomplete_left_vertex(self, X, w):
        """Given an X in L and a vertex w, returns the x in K, such that (X,x) is adjacent to w, if such x exists (In
        this case it is unique.). Otherwise returns None. (If this norm graph is bipartite, (X,x) is considered to be
        a left vertex and w is considered to be a right vertex. Note, that if this norm graph is not bipartite,
        then edgecomplete_left_vertex(X, w) == edgecomplete_right_vertex(w, X).) """
        h = self.op_norm(X, w[0])
        if self.K.is_zero(h):
            return None
        h -= w[1]
        if self.variant == NormGraphVariant.NG and w == (X, h):
            return None
        return h

    def starcomplete_right_vertex(self, A, Y):
        """Given a vertex set A and a Y in L, returns the y in K, such that all vertices in A are adjacent to (Y,y),
        if such y exists (In this case it is unique.). Otherwise returns None. (If this norm graph is bipartite,
        A is considered to be a left vertex set and (Y,y) is considered to be a right vertex. Note, that if this norm
        graph is not bipartite, then starcomplete_right_vertex(A, Y) == starcomplete_left_vertex(Y, A).) """
        return unique_element(imap(lambda v: self.edgecomplete_right_vertex(v, Y), A))

    def starcomplete_left_vertex(self, X, B):
        """Given an X in L and a vertex set B, returns the x in K, such that (X,x) is adjacent to all vertices in B,
        if such x exists (In this case it is unique.). Otherwise returns None. (If this norm graph is bipartite, (X,
        x) is considered to be a left vertex and B is considered to be a right vertex set. Note, that if this norm
        graph is not bipartite, then starcomplete_left_vertex(X, B) == starcomplete_right_vertex(B, X).) """
        return unique_element(imap(lambda v: self.edgecomplete_left_vertex(X, v), B))

    # types of pairs

    def type_equivalence_class(self, t):
        C = self.K.frobenius_equivalence_class(t)
        if -t in C:
            return C
        else:
            return C | self.K.frobenius_equivalence_class(-t)

    def type_equivalence_classes(self):
        H = set()
        R = []
        for x in self.K.unit_group:
            if x not in H:
                C = self.type_equivalence_class(x)
                H |= C
                R.append(C)
        return R

    # count subgraphs
    def count_krs(self, r, s):
        from projective_normgraph.degree_distribution.strategy import FullGenericDDCS
        degree_distribution = self.degree_distribution(r, FullGenericDDCS(), cache_load=True)
        distribution = degree_distribution.distribution()
        return sum(distribution[i] * sage.all.binomial(i, s) for i in distribution) * (1 if r == s else 2)

    def count_krs_estimating(self, r, s, sample_size=1000):
        from projective_normgraph.degree_distribution.strategy import RandomGenericSampleDDCS
        degree_distribution = self.degree_distribution(r, RandomGenericSampleDDCS(sample_size), cache_load=True)
        S = self.num_generic_vertex_sets(r) / degree_distribution.size()
        distribution = degree_distribution.distribution()
        return float(sum(S * distribution[i] * sage.all.binomial(i, s) for i in distribution) * (1 if r == s else 2))

    def count_subgraph_copies_bruteforce(self, H):
        """"Returns the number of copies of H in this graph using a brute-force method."""
        if not isinstance(H, sage.graphs.graph.Graph):
            raise ValueError("H is not a graph.")
        return self.sage_graph().subgraph_search_count(H) // H.automorphism_group().order()

    def count_krs_bruteforce(self, r, s):
        """"Returns the number of copies of K_{r,s} in this graph using a brute-force method."""
        return self.count_subgraph_copies_bruteforce(sage.all.graphs.CompleteBipartiteGraph(r, s))

    def contains_isomorphic_subgraph_bruteforce(self, H):
        """"Decides whether this graph contains a copy of H using a brute-force method."""
        if not isinstance(H, sage.graphs.graph.Graph):
            raise ValueError("H is not a graph.")
        try:
            self.sage_graph().subgraph_search_iterator(H)
            return True
        except StopIteration:
            return False


class BipartiteNormGraph(NormGraph):

    def __init__(self, q, t):
        super(BipartiteNormGraph, self).__init__(q, t, NormGraphVariant.BMNG)

    @staticmethod
    def it(t, min_q=3, max_q=sage.all.infinity, primes_only=False):
        return NormGraph.it(t, NormGraphVariant.BMNG, min_q, max_q, primes_only)
