# -*- coding: utf-8 -*-
import logging
import os
from itertools import cycle

import sage.all

from projective_normgraph import config

logger = logging.getLogger(__name__)


class MinusInfinity(object):

    def __repr__(self):
        return "-∞"

    def __eq__(self, other):
        return isinstance(other, MinusInfinity)


class ZechLogarithmCache(object):

    def __init__(self, finite_field):
        self.finite_field = finite_field
        self._cache_file = "ZechLogarithm_{}_{}.sobj".format(self.finite_field.p, self.finite_field.k)

        self.__load_cache()

    def __get_cache_file(self):
        cache_dir = config.get_cache_dir()
        if not os.path.isdir(cache_dir):
            raise ValueError("Cache directory does not exist.")
        return os.path.join(cache_dir, self._cache_file)

    def __load_cache(self):
        cache_file = self.__get_cache_file()
        if os.path.isfile(cache_file):
            logger.debug("Loading Zech logarithm cache for {!r}".format(self.finite_field))
            self._cache = sage.all.load(cache_file)
            logger.debug("Cache loaded")
        else:
            logger.debug("Building Zech logarithm cache...")
            self._cache = self.__build_cache()
            logger.debug("Saving cache...")
            sage.all.save(self._cache, self.__get_cache_file())
            logger.debug("Cache saved.")

    def __build_cache(self):
        # Build cache for field
        log_cache = {self.finite_field.field(0): self.finite_field.MI}
        h = self.finite_field.field(1)
        for i in self.finite_field.unit_group:
            log_cache[h] = i
            h *= self.finite_field.field_generator

        # Build Zech cache
        zech_cache = [None] * (self.finite_field.q - 1)
        h = self.finite_field.field(1)
        for i in xrange(self.finite_field.q - 1):
            zech_cache[i] = log_cache[h + 1]
            h *= self.finite_field.field_generator

        return zech_cache

    def __getitem__(self, key):
        return self._cache[key]


class FiniteField(object):

    @staticmethod
    def of_size(q):
        factors = sage.all.Integer(q).factor()
        if len(factors) != 1:
            raise ValueError("q={} is not a prime power.".format(q))
        return FiniteField(factors[0][0], factors[0][1])

    def __repr__(self):
        """Returns the symbolic name of this finite field"""
        if self.k == 1:
            return "FiniteField({})".format(self.q)
        else:
            return "FiniteField({}={}^{})".format(self.q, self.p, self.k)

    def __init__(self, p, k):
        if not sage.all.is_prime(p):
            raise ValueError("p must be prime.")
        if k <= 0:
            raise ValueError("The exponent must be positive.")
        self.p = p
        self.k = k
        self.q = p ** k

        # The MinusInfinity-Singleton
        self.MI = MinusInfinity()

        self.field = sage.all.GF(self.q, 'a')
        self.field_generator = self.field.multiplicative_generator()

        self.unit_group = sage.all.Integers(self.q - 1)
        self.e = 0 if self.p == 2 else self.unit_group((self.q - 1) // 2)

        self.cache = None

    def load_cache(self):
        if self.cache is None:
            self.cache = ZechLogarithmCache(self)

    def size(self):
        return self.q

    def characteristic(self):
        return self.p

    def iterator(self):
        """Returns a generator yielding all elements of the large field."""
        yield self.MI
        for a in self.unit_group:
            yield a

    def zero(self):
        return self.MI

    def is_zero(self, element):
        return element == self.MI

    def one(self):
        return self.unit_group.zero()

    def is_one(self, element):
        return element == self.unit_group.zero()

    def minus_one(self):
        return self.one() if self.p == 2 else self.e

    def generator(self):
        return self.unit_group.one()

    def random_element(self):
        h = sage.all.ZZ.random_element(0, self.q)
        return self.zero() if h == self.q - 1 else self.unit_group(h)

    def zech_logarithm(self, a):
        return self.MI if a == self.e else self.cache[a]

    def sum(self, a, b):
        if a == self.MI:
            return b
        elif b == self.MI:
            return a
        else:
            d = b - a
            return self.MI if d == self.e else a + self.cache[d]

    def neg(self, a):
        if a == self.MI:
            return a
        else:
            return a + self.e

    def product(self, a, b):
        if a == self.MI or b == self.MI:
            return self.MI
        else:
            return a + b

    def inv(self, a):
        if a == self.MI:
            raise ValueError("Zero has no multiplicative inverse.")
        return -1 * a

    def difference(self, a, b):
        if b == self.MI:
            return a
        else:
            return self.sum(a, b + self.e)

    def quotient(self, a, b):
        return self.product(a, self.inv(b))

    def power(self, a, k):
        if a == self.MI:
            return self.MI
        else:
            return k * a

    def frobenius(self, a):
        return self.power(a, self.p)

    def frobenius_equivalence_class(self, a):
        equivalence_class = {a}
        h = self.frobenius(a)
        while h != a:
            equivalence_class.add(a)
            a = self.frobenius(a)
        return equivalence_class

    def __eq__(self, other):
        return isinstance(other, FiniteField) and other.q == self.q

    def __len__(self):
        return self.q

    def __iter__(self):
        return self.iterator()


class ZechLogarithmNormCache(object):

    def __init__(self, finite_field_extension):
        self.ffe = finite_field_extension
        self._cache_file = "ZechLogarithmNorm_{}_{}_{}.sobj".format(
            finite_field_extension.L.p,
            finite_field_extension.L.k,
            finite_field_extension.K.k
        )

        self.__load_cache()

    def __get_cache_file(self):
        cache_dir = config.get_cache_dir()
        if not os.path.isdir(cache_dir):
            raise ValueError("Cache directory does not exist.")
        return os.path.join(cache_dir, self._cache_file)

    def __load_cache(self):
        cache_file = self.__get_cache_file()
        if os.path.isfile(cache_file):
            logger.debug("Loading Zech logarithm norm cache for {!r}".format(self.ffe))
            self._cache = sage.all.load(cache_file)
            logger.debug("Cache loaded")
        else:
            logger.debug("Building Zech logarithm norm cache...")
            self._cache = self.__build_cache()
            logger.debug("Saving cache...")
            sage.all.save(self._cache, self.__get_cache_file())
            logger.debug("Cache saved.")

    def __build_cache(self):
        # Build cache for field
        log_cache = {self.ffe.L.field(0): self.ffe.K.MI}
        h = self.ffe.L.field(1)
        for i in cycle(self.ffe.K.unit_group):
            log_cache[h] = i
            h *= self.ffe.L.field_generator
            if h == self.ffe.L.field(1):
                break

        # Build Zech norm cache
        zech_norm_cache = [None] * (self.ffe.L.q - 1)
        h = self.ffe.L.field(1)
        for i in xrange(self.ffe.L.q - 1):
            zech_norm_cache[i] = log_cache[h + 1]
            h *= self.ffe.L.field_generator

        return zech_norm_cache

    def __getitem__(self, key):
        return self._cache[key]


class FiniteFieldExtension(object):

    @staticmethod
    def of_sizes(q, p):
        q_factors = sage.all.Integer(q).factor()
        if len(q_factors) != 1:
            raise ValueError("q={} is not a prime power.".format(q))
        p_factors = sage.all.Integer(p).factor()
        if len(p_factors) != 1:
            raise ValueError("p={} is not a prime power.".format(p))
        if q % p != 0:
            raise ValueError("p={} does not divide q={}.".format(p, q))
        return FiniteFieldExtension(q_factors[0][0], q_factors[0][1], p_factors[0][1])

    def __repr__(self):
        """Returns the symbolic name of this finite field extension"""
        return "FiniteFieldExtension({} over {})".format(self.L, self.K)

    def __init__(self, p, k, m):
        if not sage.all.is_prime(p):
            raise ValueError("p must be prime.")
        if k <= 0 or m <= 0:
            raise ValueError("The exponents must be positive.")
        if k % m != 0:
            raise ValueError("The first exponent must be divisible by the second.")
        self.L = FiniteField(p, k)
        self.K = FiniteField(p, m)

        self.cache = None

    def load_cache(self):
        if self.cache is None:
            self.cache = ZechLogarithmNormCache(self)

    def degree(self):
        return self.L.k // self.K.k

    def norm(self, a):
        if a == self.L.MI:
            return self.K.MI
        else:
            return self.K.unit_group(a)

    def embedding(self, a):
        if a == self.K.MI:
            return self.L.MI
        else:
            return self.L.unit_group(a) * ((self.L.q - 1) // (self.K.q - 1))

    def zech_logarithm_norm(self, a):
        return self.K.MI if a == self.L.e else self.cache[a]

    def sum_norm(self, a, b):
        if a == self.L.MI:
            return self.norm(b)
        elif b == self.L.MI:
            return self.norm(a)
        else:
            d = b - a
            return self.K.MI if d == self.L.e else self.K.unit_group(a) + self.cache[d]

    def difference_norm(self, a, b):
        if b == self.L.MI:
            return a
        else:
            return self.sum_norm(a, b + self.L.e)

    def norm_preimage(self, a):
        if a == self.K.MI:
            yield self.L.MI
        else:
            for x in xrange(a, self.L.q, self.K.q - 1):
                yield self.L.unit_group(x)

    def __eq__(self, other):
        return isinstance(other, FiniteFieldExtension) and other.K == self.K and other.L == self.L
