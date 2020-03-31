# -*- coding: utf-8 -*-

import sage.all

import os
import re
from time import strftime


class DegreeDistribution(object):

    def __init__(self, ng, r, sample_description):
        self.ng = ng
        self.r = r
        self.sample_description = sample_description
        self.distribution = {}

    def __str__(self):
        """Returns a description of this degree distribution."""
        return "Distribution of degrees of {} {}-vx-sets ({}) in {!r}".format(self.size(), self.r,
                                                                              self.sample_description, self.ng)

    def __repr__(self):
        return "DegreeDistribution({}, {}, {})".format(self.ng, self.r, self.sample_description)

    def __getstate__(self):
        return self.ng, self.r, self.sample_description, self.distribution

    def __setstate__(self, state):
        self.ng, self.r, self.sample_description, self.distribution = state

    def size(self):
        return sum(self.distribution.values())

    def distribution(self):
        return self.distribution

    def inc(self, i, n=1):
        if i in self.distribution:
            self.distribution[i] += n
        else:
            self.distribution[i] = n

    def mul(self, n):
        for i in self.distribution:
            self.distribution[i] *= n

    def merge(self, degree_distribution):
        other = degree_distribution.distribution
        for i in other:
            if i in self.distribution:
                self.distribution[i] += other[i]
            else:
                self.distribution[i] = other[i]

    def __eq__(self, other):
        return self.__getstate__() == other.__getstate__()

    def __ne__(self, other):
        return not self.__eq__(other)

    @staticmethod
    def load(ng, sample_description, r, results_dir):
        expr = re.compile(r'(.*){!s}\({},{}\)_{}\({}\).sobj'.format(ng.variant, ng.q, ng.t, r, sample_description))
        for file_name in os.listdir(results_dir):
            if expr.match(file_name):
                yield sage.all.load(os.path.join(results_dir, file_name))

    def save(self, sample_description, results_dir):
        filename = "{}_{!r}_{}({}).sobj".format(strftime("%Y-%m-%d_%H:%M:%S"), self.ng, sample_description, self.r)
        sage.all.save(self, os.path.join(results_dir, filename))

    def print_statistics(self):
        total = sum(self.distribution.values())
        print sage.all.table(
            [("size", "count", "share")]
            + map(lambda i: (i, self.distribution[i], self.distribution[i] / float(total)), sorted(self.distribution))
            + [("total", total, 1)],
            header_row=True, header_column=True)
