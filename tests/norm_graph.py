from projective_normgraph.normgraph import *

# Graph parameters
for variant in NormGraphVariant:
    for t in xrange(2, 7):
        if variant == NormGraphVariant.MNG and t % 2 == 0:
            continue
        for A in NormGraph.it(t, variant):
            if A.num_vertices() > 200:
                break
            print "Testing {}".format(A)
            g = A.sage_graph()
            assert g.num_verts() == A.num_vertices()
            assert g.num_edges() == A.num_edges()
            assert g.average_degree() == A.average_degree()
            assert g.density() == A.density()
            assert g.automorphism_group().order() == A.automorphism_group_order()
