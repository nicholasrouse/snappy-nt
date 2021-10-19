#  snappy-nt: extensions to SnapPy for arithmetic invariants


## Introduction

This package allows for easily computing the (invariant) trace field, (invariant) quaternion algebra, and denominators associated to finite volume hyperbolic 3-manifolds. It must be run in [SageMath](https://www.sagemath.org/) with [SnapPy](http://snappy.math.uic.edu/index.html) installed. The package depends on the excellent SnapPy package maintained by Culler, Dunfield, and Goerner but is a separate project.

The project is still in its early stages, so the documentation is incomplete, but one can find installation and basic usage instructions below.

## Installation

Get the repository with

    $ git clone https://github.com/nicholasrouse/snappy-nt

To install into SageMath, first change directories using `cd snappy-nt`. Then install the package using

    $ sage --python -m pip install -e .

It's a good idea to check that snappy is installed.

    $ sage --pip install snappy

If you run into issues with compatibility on your system, it might be a good idea to use the SnapPy maintainer's [Docker image](https://github.com/3-manifolds/sagedocker) on which I have tested the package.

## Usage

To start SageMath, usually all one needs to do is use the `sage` command. Inside SageMath, one needs to import the package each time using

    sage: from snappynt.ManifoldNT import ManifoldNT

Creating manifolds works similarly to creating them in SnapPy:

    sage: mfld = ManifoldNT('4_1')
    sage: mfld.compute_arithmetic_invariants()
    sage: mfld.p_arith()
    Orbifold name: 4_1(0,0)
    Volume: 2.02988321281931
    Trace field: Number Field in z with defining polynomial x^2 - x + 1 with z = 0.50000000000000000? + 0.866025403784439?*I
            Signature: (0, 1)
            Discriminant: -3
    Quaternion algebra:
            Hilbert symbol: (4*z - 8, -z)
            Finite ramification: set()
            Finite ramification residue characteristics: Counter()
            Number of ramified real places: 0
            Ramified real places: []
    Invariant Trace field: Number Field in z with defining polynomial x^2 - x + 1 with z = 0.50000000000000000? + 0.866025403784439?*I
            Signature: (0, 1)
            Discriminant: -3
    Invariant quaternion algebra:
            Hilbert symbol: (-32*z + 16, 12*z)
            Finite ramification: set()
            Finite ramification residue characteristics: Counter()
            Number of ramified real places: 0
            Ramified real places: []
    Integer traces: True
    Arithmetic: True

The main arithmetic invariants the package computes are the (invariant) trace field, the (invariant) quaternion algebra, and denominators. They can be computed by

    sage: mfld = ManifoldNT('m010(-1,2)')
    sage: mfld.trace_field()
    Number Field in z with defining polynomial x^4 - 2*x^2 + 2 with z = -1.098684113467810? - 0.4550898605622274?*I
    sage: mfld.invariant_trace_field()
    Number Field in z with defining polynomial x^2 + 1 with z = -1*I
    sage: mfld.quaternion_algebra()
    Quaternion Algebra (-2*z^2, -3*z^2 + 2) with base ring Number Field in z with defining polynomial x^4 - 2*x^2 + 2 with z = -1.098684113467810? - 0.4550898605622274?*I
    sage: mfld.invariant_quaternion_algebra()
    Quaternion Algebra (-8, -12*z + 4) with base ring Number Field in z with defining polynomial x^2 + 1 with z = -1*I
    sage: mfld.denominators()
    set()

One can try to compute all the arithmetic invariants at once with the `compute_arithmetic_invariants()` method and print them with the `p_arith()` method.

    sage: mfld = ManifoldNT('m003(-2,3)')
    sage: mfld.compute_arithmetic_invariants()
    sage: mfld.p_arith()
    Orbifold name: m003(-2,3)
    Volume: 0.981368828892232
    Trace field: Number Field in z with defining polynomial x^4 - x - 1 with z = -0.2481260628026220? + 1.033982060975968?*I
            Signature: (2, 1)
            Discriminant: -283
    Quaternion algebra:
            Hilbert symbol: (-z^3 + z^2 - 3, -z^3 - z^2 + z)
            Finite ramification: set()
            Finite ramification residue characteristics: Counter()
            Number of ramified real places: 2
            Ramified real places: ['z |--> 1.220744084605759475361685349109']
    Invariant Trace field: Number Field in z with defining polynomial x^4 - x - 1 with z = -0.2481260628026220? + 1.033982060975968?*I
            Signature: (2, 1)
            Discriminant: -283
    Invariant quaternion algebra:
            Hilbert symbol: (3*z^3 - 3*z^2 - z - 2, 4*z^3 + z^2 - 5*z - 3)
            Finite ramification: set()
            Finite ramification residue characteristics: Counter()
            Number of ramified real places: 2
            Ramified real places: ['z |--> 1.220744084605759475361685349109']
    Integer traces: True
    Arithmetic: True

For more complicated manifolds, the first pass of computing an invariant might not find it because the precision or degree used for the [LLL algorithm](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm) is too small. One can pass in explicit precision and degree choices as keyword arguments.

    sage: mfld = ManifoldNT('4_1(15,7)')
    sage: mfld.trace_field() is None
    True # The field is not found.
    sage: m.trace_field(degree=26, prec=10000)
    Number Field in z with defining polynomial x^26 - 2*x^25 - 5*x^24 + 15*x^23 + 15*x^22 - 72*x^21 - 51*x^20 + 309*x^19 + 121*x^18 - 1091*x^17 - 82*x^16 + 3049*x^15 - 530*x^14 - 6675*x^13 + 2533*x^12 + 10642*x^11 - 5799*x^10 - 10963*x^9 + 7287*x^8 + 6273*x^7 - 4655*x^6 - 1595*x^5 + 1249*x^4 + 135*x^3 - 111*x^2 + 5*x + 1 with z = 1.009043443852433? + 1.675767880830078?*I

Alternatively, one can just retry computing the invariant without explicit keyword arguments, and the package will attempt to guess appropriate precision and degree.

    sage: mfld = ManifoldNT('4_1(15,7)')
    sage: mfld.trace_field() # No output
    sage: mfld.trace_field() # Again no output
    sage: mfld.trace_field() # Third time's the charm
    Number Field in z with defining polynomial x^26 - 2*x^25 - 5*x^24 + 15*x^23 + 15*x^22 - 72*x^21 - 51*x^20 + 309*x^19 + 121*x^18 - 1091*x^17 - 82*x^16 + 3049*x^15 - 530*x^14 - 6675*x^13 + 2533*x^12 + 10642*x^11 - 5799*x^10 - 10963*x^9 + 7287*x^8 + 6273*x^7 - 4655*x^6 - 1595*x^5 + 1249*x^4 + 135*x^3 - 111*x^2 + 5*x + 1 with z = 1.009043443852433? + 1.675767880830078?*I

## Contact

If you would like to reach me about this project, you can open an issue, or send an email to `snappynumbertheory@outlook.com`.

## Known issues

If one tries to run

    $ sage --python
    >>> import snappynt.QuaternionAlgebraNF

one will get a mysterious import error. The solution is to either run sage rather than
sage's python, or first import sage.all. See https://groups.google.com/g/sage-support/c/3lXInDgHJe4/m/SdOpuuOkqPUJ
for more details.
