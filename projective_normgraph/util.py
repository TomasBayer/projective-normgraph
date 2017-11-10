# -*- coding: utf-8 -*-
import itertools


def all_distinct(iterable):
    """Decides whether all elements in the iterable are distinct."""
    p = set()
    for x in iterable:
        if x in p:
            return False
        p.add(x)
    return True


def unique_element(iterable):
    """If all elements in the iterable are equal, returns the unique element, otherwise returns None."""
    x = next(iterable)
    return x if all(x == y for y in iterable) else None


def produce(function, n=None):
    """Repeatedly invokes the function without any arguments and yields its return value. If n is given,
    stops after n iterations."""
    if n is None:
        while True:
            yield function()
    else:
        for _ in xrange(n):
            yield function()


def truncate(iterable, n):
    """Returns an iterator which yields the first n elements of the given iterable."""
    return itertools.islice(iterable, 0, n)


def chunk(iterable, chunk_size):
    it = iter(iterable)
    return iter(lambda: tuple(itertools.islice(it, chunk_size)), ())


def skip_until(iterable, o):
    it = iter(iterable)
    for x in it:
        if x == o:
            break
    for x in it:
        yield x


def prime_powers(start, stop=None):
    """Functional analogue to sage.rings.arith.primes() for prime powers."""
    from sage.rings.infinity import infinity
    from sage.rings.integer_ring import ZZ

    start = ZZ(start)
    if stop is None:
        stop = start
        start = ZZ(2)
    elif stop != infinity:
        stop = ZZ(stop)
    n = start - 1
    while True:
        n = n.next_prime_power()
        if n < stop:
            yield n
        else:
            return
