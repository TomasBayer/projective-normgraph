from projective_normgraph.finitefield import *

to_field = lambda F, a: F.field(0) if F.is_zero(a) else F.field_generator ** a

TEST_LENGTH = 100

for p in sage.all.primes_first_n(10):
    for k in xrange(1, 5):
        F = FiniteField(p, k)
        print "Testing {}".format(F)
        F.load_cache()

        to_F = lambda a: to_field(F, a)

        for _ in xrange(TEST_LENGTH):
            a = F.random_element()
            b = F.random_element()
            assert to_F(F.sum(a, b)) == to_F(a) + to_F(b), "sum(a, b) wrong for a={}, b={}".format(a, b)
            assert to_F(F.neg(a)) == -1 * to_F(a), "neg(a) wrong for a={}".format(a)
            assert to_F(F.difference(a, b)) == to_F(a) - to_F(b), "difference(a, b) wrong for a={}, b={}".format(a, b)
            assert to_F(F.product(a, b)) == to_F(a) * to_F(b), "product(a, b) wrong for a={}, b={}".format(a, b)
            if b != F.MI:
                assert to_F(F.inv(b)) == 1 / to_F(b), "inv(a) wrong for a={}".format(a)
                assert to_F(F.quotient(a, b)) == to_F(a) / to_F(b), "quotient(a, b) wrong for a={}, b={}".format(a, b)

        for l in sage.all.divisors(k):
            FE = FiniteFieldExtension(p, k, l)
            print "Testing {}".format(FE)

            FE.load_cache()
            FE.L.load_cache()
            FE.K.load_cache()

            to_L = lambda a: to_field(FE.L, a)
            to_K = lambda a: to_field(FE.K, a)
            embedding = FE.L.field.subfields(l, 'a')[0][1]
            norm = lambda X: X ** ((p ** k - 1) // (p ** l - 1))

            # By SageMath specifications, GF is implemented using Conway polynomials
            # and the multiplicative_generator()-method returns the lexicographically
            # smallest generator. It follows from the definition that the field norm
            # applied on the generators of the large field as provided by that
            # method always yields the generator of the smaller field, which is provided
            # by the same method applied on the smaller field.
            #
            # This fact is used by the Zech logarithm cache builder.
            #
            # Sanity check:
            assert embedding(FE.K.field_generator) == norm(
                FE.L.field_generator), "The generators are not norm-compatible!"

            for _ in xrange(TEST_LENGTH):
                a = F.random_element()
                assert embedding(to_K(FE.norm(a))) == norm(to_L(a)), "norm(a) wrong for a={}".format(a)

            for _ in xrange(TEST_LENGTH):
                a = F.random_element()
                assert embedding(to_K(a)) == to_L(FE.embedding(a)), "embedding(a) wrong for a={}".format(a)

            for _ in xrange(TEST_LENGTH):
                a = F.random_element()
                b = F.random_element()
                assert FE.sum_norm(a, b) == FE.norm(FE.L.sum(a, b)), "sum_norm(a, b) wrong for a={}, b={}".format(a, b)

            for _ in xrange(TEST_LENGTH):
                a = F.random_element()
                b = F.random_element()
                assert FE.difference_norm(a, b) == FE.norm(FE.L.difference(a, b)), \
                    "difference_norm(a, b) wrong for a={}, b={}".format(a, b)
