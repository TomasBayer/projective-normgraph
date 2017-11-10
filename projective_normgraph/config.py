# -*- coding: utf-8 -*-
import multiprocessing
import os

CPU_COUNT = multiprocessing.cpu_count()


class Config(object):

    def __init__(self):
        self._cache_dir = "cache/"
        self._results_dir = "results/"
        self._max_processes = 1

    def get_cache_dir(self):
        return self._cache_dir

    def set_cache_dir(self, value):
        cache_dir = os.path.expanduser(value)
        if not os.path.isdir(cache_dir):
            raise ValueError("cache_dir does not exist.")
        self._cache_dir = cache_dir

    def get_results_dir(self):
        return self._results_dir

    def set_results_dir(self, value):
        results_dir = os.path.expanduser(value)
        if not os.path.isdir(results_dir):
            raise ValueError("results_dir does not exist.")
        self._results_dir = results_dir

    def get_max_processes(self):
        return self._max_processes

    def set_max_processes(self, value):
        if value < 1:
            raise ValueError("max_processes has to be 1 at least.")
        if value > CPU_COUNT:
            raise ValueError("max_processes should not be larger than the number of CPUs ({}).".format(CPU_COUNT))
        self._max_processes = value


config = Config()
