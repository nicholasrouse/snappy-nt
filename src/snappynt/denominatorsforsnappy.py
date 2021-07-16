"""
Ideally this module will get incorporated into other at some point. It was written
before ManifoldNT was, so we worked systematically with SnapPy's Manifold objects
rather than a custom subclass.
"""

# For some reason I thought this was wrong at some point.
# Perhaps because I thought the generators corresponded to matrix entries.
# According to SnapPy source code and doc strings, trace_field_gens()
# Should give generators for the trace field.
# Need to make sure this is reliable for (non)integral traces.
# NMR Jul-05-2020

import functools
import itertools
import math
import shelve
import time

import snappy
from sage.all import denominator, factor

# import large_knots


def find_prime_factors_in_a_set(aSet):  # aSet should be (positive) integers
    primes = set()
    for element in aSet:
        factorization = list(factor(element))
        for (prime, power) in factorization:
            primes.add(prime)
    return primes


def find_denominators_for_a_manifold(
    manifold, prec=10 ** 4, degree=10, verbosity=False
):
    unfactored_denoms = set()
    field_generators = manifold.trace_field_gens().find_field(
        prec=prec, degree=degree, optimize=True
    )
    if not field_generators:
        return None
        if verbosity:
            print(f"Unable to find trace field at precision={prec} and degree={degree}")
    for element in field_generators[2]:
        unfactored_denoms.add(denominator(element.norm()))
    primes = find_prime_factors_in_a_set(unfactored_denoms)
    return primes


def find_denominators_for_a_manifold_variable_precision(
    manifold,
    starting_prec=10 ** 3,
    starting_degree=10,
    prec_increment=10 ** 3,
    degree_increment=5,
    max_prec=10 ** 5,
    max_degree=100,
    trace=False,
):
    prec = starting_prec
    degree = starting_degree
    primes = None
    while primes is None:
        if trace:
            print(
                str(manifold) + ":", f"Trying with precision={prec} and degree={degree}"
            )
        primes = find_denominators_for_a_manifold(manifold, prec, degree)
        if prec == max_prec and degree == max_degree:
            return primes
        if prec + prec_increment <= max_prec:
            prec = prec + prec_increment
        else:
            prec = max_prec
        if degree + degree_increment <= max_degree:
            degree = degree + degree_increment
        else:
            degree = max_degree
    return primes


def pretty_seconds(seconds):  # Seems like there should be a better way?
    years = math.floor(seconds / (60 * 60 * 24 * 365))
    seconds = seconds - years * 60 * 60 * 24 * 365
    days = math.floor(seconds / (60 * 60 * 24))
    seconds = seconds - days * 60 * 60 * 24
    hours = math.floor(seconds / (60 * 60))
    seconds = seconds - hours * 60 * 60
    minutes = math.floor(seconds / (60))
    seconds = seconds - minutes * 60
    if years:
        return f"{years} years, {days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
    if days:
        return f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
    if hours:
        return f"{hours} hours, {minutes} minutes, {seconds} seconds"
    if minutes:
        return f"{minutes} minutes, {seconds} seconds"
    return f"{seconds} seconds"


# This should maybe be in a separate module. Jul-19-20.
"""
def trawl_census_iterator(
    iterator, function, start=0, stop=None, output_file=None, **kwargs
):
    if output_file is None:
        iterator_name = str(iterator)
        function_name = function.__name__
        options_string = ""
        for key in kwargs:
            key_string = str(key)
            value_string = str(value)
            options_string = options_string + f"{key_string}={value_string},"
        output_file = f"{iterator_name}{function_name}{options_string}"

    sliced_iterator = itertools.islice(iterator, start, stop, 1)
    sliced_iterator = list(
        sliced_iterator
    )  # This is not actually a generator. It's just a trick to get the list to be garbage collected.
    iterator_size = len(sliced_iterator)
    sliced_iterator = iter(
        sliced_iterator
    )  # Resets the iterator. Should garbage collect the list.
    running_timer = 0  # This is a running timer of how many have been computed so far.
    elements_considered = 0
    print("Finished startup tasks...")
    print("Going to consider", iterator_size, "elements")
    with shelve.open(output_file) as dictionary_name:
        for element in sliced_iterator:
            before_function_timer = time.monotonic()
            print("Going to consider", str(element))
            try:
                function_output = function(element)
            except ValueError as error:
                if not element.verify_hyperbolicity():
                    print("Considered manifold is not hyperbolic.")
                    function_output = "None"
                else:
                    print("Snappy ran into a problem with", element)
                    print(error)
                    function_output = None
            elements_considered += 1
            dictionary_name[str(element)] = function_output
            print(element, "evaluated to", str(function_output))
            print(f"Considered {elements_considered} so far.")
            after_function_timer = time.monotonic()
            one_operation_time = after_function_timer - before_function_timer
            running_timer += one_operation_time
            print(
                "It took",
                pretty_seconds(one_operation_time),
                "to compute",
                str(element),
            )
            print("Time elapsed so far:", pretty_seconds(running_timer))
            estimated_time_remaining = (running_timer / elements_considered) * (
                iterator_size - elements_considered
            )
            print("Estimated time remaining:", pretty_seconds(estimated_time_remaining))
"""

"""
def find_denominators_for_large_knots():
    i = 0
    with shelve.open("Large_Nonintegral_Knots") as nonintegraldict:
        G = large_knots.large_knots()
        for knot in G:
            datadict = dict()
            try:
                if i % 10 == 0:
                    print(i)
                i += 1
                denoms = find_denominators_for_a_manifold_variable_precision(
                    knot,
                    starting_prec=1000,
                    starting_degree=10,
                    max_prec=10 ** 5,
                    max_degree=100,
                )
                if denoms:
                    print(knot, "is nonintegral at", denoms)
                    datadict["denominators"] = denoms
                    # datadict['field info'] = manifold.trace_field_gens().find_field(prec, degree, optimize=True)
                    nonintegraldict[str(knot)] = datadict
                if denoms is None:
                    print("Couldn't find trace field for", knot)
            except ValueError as error:
                print("Snappy ran into a problem with", knot)
                print(error)
"""


def test():
    i = 0
    with shelve.open("10Crossing2ComponentLinkNonintegralData") as nonintegraldict:
        for link in snappy.LinkExteriors(crossings=10, num_cusps=2):
            datadict = dict()
            try:
                if i % 10 == 0:
                    print(i)
                i += 1
                denoms = find_denominators_for_a_manifold_variable_precision(
                    link,
                    starting_prec=1000,
                    starting_degree=10,
                    max_prec=10 ** 5,
                    max_degree=100,
                )
                if denoms:
                    print(link, "is nonintegral at", denoms)
                    datadict["denominators"] = denoms
                    templink = snappy.Link(str(link)[:-10])
                    datadict["linking_number"] = templink.linking_number()
                    datadict["DT_code"] = link.DT_code(alpha=True)
                    # datadict['field info'] = manifold.trace_field_gens().find_field(prec, degree, optimize=True)
                    nonintegraldict[str(link.identify()[1])] = datadict
                if denoms is None:
                    print("Couldn't find trace field for", link)
            except ValueError as error:
                print("Snappy ran into a problem with", link)
                print(error)


def test2(crossings, start=0, stop=None):
    i = 0
    start_timer = time.monotonic()
    running_total = 0
    dictionaryname = str(crossings) + "Crossing2ComponentLinkNonintegralData"
    linkiterator = itertools.islice(
        snappy.HTLinkExteriors(crossings=crossings, num_cusps=2), start, stop, 1
    )
    iteratorsize = len(list(linkiterator))
    linkiterator = itertools.islice(
        snappy.HTLinkExteriors(crossings=crossings, num_cusps=2), start, stop, 1
    )  # Resets the iterator.
    print(iteratorsize)
    candidatesconsidered = 0
    with shelve.open(dictionaryname) as nonintegraldict:
        for link in linkiterator:
            if candidatesconsidered % 100 == 0:
                print("Considered", candidatesconsidered, "links so far.")
            datadict = dict()
            G = link.fundamental_group()
            N = len(G.generators())
            if N == 2:
                try:
                    if i != 0:
                        newtimer = time.monotonic()
                        last_timer = newtimer - start_timer
                        running_total += last_timer
                        start_timer = time.monotonic()
                        estimated_time_remaining = (iteratorsize - i) * (
                            running_total / i
                        )
                        print(
                            f"Computed {i} links so far. \nIt took "
                            + pretty_seconds(last_timer)
                            + " to compute the previous link, and "
                        )
                        print(
                            pretty_seconds(running_total)
                            + " to compute all links so far."
                        )
                        print(
                            "Estimated time remaining: "
                            + pretty_seconds(estimated_time_remaining)
                        )
                    i += 1
                    denoms = find_denominators_for_a_manifold_variable_precision(
                        link,
                        starting_prec=1000,
                        starting_degree=10,
                        max_prec=10 ** 5,
                        max_degree=100,
                    )
                    if denoms:
                        print("!!!", link, "is nonintegral at", denoms)
                        datadict["denominators"] = denoms
                        templink = snappy.Link(str(link)[:-10])
                        datadict["linking_number"] = templink.linking_number()
                        datadict["DT_code"] = link.DT_code(alpha=True)
                        # datadict['field info'] = manifold.trace_field_gens().find_field(prec, degree, optimize=True)
                        nonintegraldict[str(link.identify())] = datadict
                    if denoms is None:
                        print("Couldn't find trace field for", link)
                except ValueError as error:
                    if not link.verify_hyperbolicity():
                        print("Link is not hyperbolic.")
                    else:
                        print("Snappy ran into a problem with", link)
                        print(error)
            candidatesconsidered += 1


def find_denoms_partial(
    starting_prec=2 * 10 ** 3,
    starting_degree=10,
    max_prec=10 ** 5,
    max_degree=100,
    prec_increment=5000,
    degree_increment=5,
    trace=False,
):
    """
    Creates the relevant partial function to pass in a Manifold to look for denominators
    """
    f = functools.partial(
        find_denominators_for_a_manifold_variable_precision,
        starting_prec=starting_prec,
        starting_degree=starting_degree,
        degree_increment=degree_increment,
        max_prec=max_prec,
        prec_increment=prec_increment,
        max_degree=max_degree,
        trace=trace,
    )
    return f
