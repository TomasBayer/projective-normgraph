# -*- coding: utf-8 -*-
from multiprocessing import current_process
from sage.parallel.multiprocessing_sage import parallel_iter


class ProcessingPool(object):

    def __init__(self, max_processes):
        self.max_processes = max_processes

    def run(self, method, args_iterator):
        if self.max_processes == 1 or current_process().name != 'MainProcess':
            for args in args_iterator:
                yield (args, method(self, *args))
        else:
            for ((args, kwargs), ret) in parallel_iter(self.max_processes, method,
                                                       (((self,) + args, {}) for args in args_iterator)):
                yield (args, ret)
