# -*- coding: utf-8 -*-

import sage.all

class RandomGraph(object):
    """Probability space on the set of simple graphs with fixed vertex number and edge probability"""

    def __init__(self, num_vertices, edge_probability):
        if num_vertices < 1:
            raise ValueError("The number of vertices must be positive.")
        if edge_probability < 0:
            raise ValueError("The edge probability must not be negative.")
        if edge_probability > 1:
            raise ValueError("The edge probability must not be larger than 1.")
        self.num_vertices = num_vertices
        self.edge_probability = edge_probability

    @staticmethod
    def of(graph):
        """Returns the associated random graph.

        Args:
            graph: a sage graph object

        Returns:
            a RandomGraph(n,p) where n is the number of vertices of G and p is the edge density of G,
            i.e. e(G)/(n choose 2).
        """
        if not isinstance(graph, sage.graphs.graph.Graph):
            raise ValueError("G is not a graph.")
        return RandomGraph(graph.num_verts(), graph.num_edges() / sage.all.binomial(graph.vertices, 2))

    def __repr__(self):
        """Returns the symbolic name of this random graph"""
        return "RandomGraph({},{})".format(self.num_vertices, self.edge_probability)

    def sample(self):
        """Returns a graph randomly chosen from this probability space."""
        return sage.all.graphs.RandomGNP(self.num_vertices, self.edge_probability)

    def evaluate(self, graph_function, sample_size=100):
        """Evaluates a given function on n randomly chosen graphs and returns the average.

        Randomly chooses a graphs from this probability space n times, applies func on each of them and returns
        the average of the returned values.

        Args:
           graph_function: A function which takes exactly one argument - a graph - and returns a real number.
           sample_size: The size of the sample.

        Returns:
            The average of the values return from the function evaluations.
        """
        if sample_size < 1:
            raise ValueError("The sample size must be positive")
        results = [graph_function(self.sample()) for _ in xrange(sample_size)]
        return sum(results) / float(len(results))

    def num_vertices(self):
        """Returns the number of vertices."""
        return self.num_vertices

    def num_edges(self):
        """Returns the expected number of edges."""
        return self.edge_probability * sage.all.binomial(self.num_vertices, 2)

    def num_subgraphs_by_parameters(self, num_vertices, num_edges, num_automorphisms):
        """Returns the expected number of copies of a graph with the given parameters"""
        return sage.all.binomial(self.num_vertices, num_vertices) \
               * sage.all.factorial(num_vertices) \
               / num_automorphisms \
               * self.edge_probability ** num_edges

    def num_subgraphs(self, subgraph):
        """Returns the expected number of copies of the given graph."""
        return self.num_subgraphs_by_parameters(
            subgraph.num_verts(),
            subgraph.num_edges(),
            subgraph.automorphism_group().order()
        )

    def num_triangles(self):
        """Returns the expected number of triangles."""
        return self.num_subgraphs_by_parameters(3, 3, 6)


class RandomBipartiteGraph(object):

    def __init__(self, num_left_vertices, num_right_vertices, edge_probability):
        if num_left_vertices < 1:
            raise ValueError("The number of left vertices must be positive.")
        if num_right_vertices < 1:
            raise ValueError("The number of right vertices must be positive.")
        if edge_probability < 0:
            raise ValueError("The edge probability must not be negative.")
        if edge_probability > 1:
            raise ValueError("The edge probability must not be larger than 1.")
        self.num_left_vertices = num_left_vertices
        self.num_right_vertices = num_right_vertices
        self.edge_probability = edge_probability

    @staticmethod
    def of(graph):
        """Returns the associated bipartite random graph.

        Args:
            graph: a sage graph object containing a bipartite graph

        Returns:
             a RandomBipartiteGraph(l,r,p) where l and r are the number of vertices in the left and right part of the
             vertex partition of G respectively and p is the proportion of edges in between the parts, i.e. e(G)/(l*r).
        """
        if not isinstance(graph, sage.graphs.graph.Graph):
            raise ValueError("G is not a graph.")
        if not graph.is_bipartite():
            raise ValueError("G is not a bipartite graph.")
        bipartite_graph = sage.all.graphs.BipartiteGraph(graph)
        num_left_vertices, num_right_vertices = len(bipartite_graph.left), len(bipartite_graph.right)
        return RandomBipartiteGraph(
            num_left_vertices,
            num_right_vertices,
            bipartite_graph.num_edges() / (num_left_vertices * num_right_vertices)
        )

    def __repr__(self):
        """Returns the symbolic name of this random graph"""
        return "RandomBipartiteGraph({},{},{})".format(self.num_left_vertices, self.num_right_vertices,
                                                       self.edge_probability)

    def sample(self):
        """Returns a graph randomly chosen from this probability space."""
        return sage.all.graphs.BipartiteGraph(
            sage.all.graphs.RandomBipartite(self.num_left_vertices, self.num_right_vertices, self.edge_probability))

    def evaluate(self, graph_function, sample_size=100):
        """Evaluates a given function on n randomly chosen graphs and returns the average.

        Randomly chooses n graphs from this probability space, applies func on each of them and return the average of
        the result values.

        Args:
           graph_function: A function which takes exactly one argument - a graph - and returns a real number.
           sample_size: The size of the sample.

        Returns:
            The average of the values return from the function evaluations.
        """
        if sample_size < 1:
            raise ValueError("The sample size must be positive")
        results = [graph_function(self.sample()) for _ in xrange(sample_size)]
        return sum(results) / float(len(results))

    def num_vertices(self):
        """Returns the number of vertices."""
        return self.num_left_vertices + self.num_right_vertices

    def num_edges(self):
        """Returns the expected number of edges."""
        return self.edge_probability * self.num_left_vertices * self.num_right_vertices

    def num_subgraphs_by_parameters(self, num_left_vertices, num_right_vertices, num_edges, num_automorphisms):
        """Returns the expected number of copies of a bipartite graph with the given parameters"""
        return sage.all.binomial(self.num_left_vertices, num_left_vertices) \
               * sage.all.factorial(num_left_vertices) \
               * sage.all.binomial(self.num_right_vertices, num_right_vertices) \
               * sage.all.factorial(num_right_vertices) \
               / num_automorphisms \
               * self.edge_probability ** num_edges

    def num_subgraphs(self, subgraph):
        """Returns the expected number of copies of the given graph."""
        if subgraph.is_bipartite():
            bipartite_graph = sage.all.graphs.BipartiteGraph(subgraph)
            return self.num_subgraphs_by_parameters(
                len(bipartite_graph.left),
                len(bipartite_graph.right),
                subgraph.num_edges(),
                subgraph.automorphism_group().order()
            )
        else:
            return 0

    def num_cycles(self, n):
        """Returns the expected number of cycles of length n."""
        if n < 3:
            raise ValueError("The length of the cycle must not be smaller than 3.")
        if n % 2 == 0:
            return self.num_subgraphs_by_parameters(n // 2, n // 2, n, 2 * n)
        else:
            return 0
