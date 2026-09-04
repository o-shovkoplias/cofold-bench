"""Unit tests for the pDockQ implementation on synthetic data."""

import math

import numpy as np

from cofoldbench.pdockq import DEFAULT, interface_stats, pdockq_from_arrays, pdockq_from_x


def test_sigmoid_limits():
    # far below x0 -> b ; far above -> L + b
    assert math.isclose(pdockq_from_x(-1e6), DEFAULT.b, rel_tol=1e-9)
    assert math.isclose(pdockq_from_x(1e6), DEFAULT.L + DEFAULT.b, rel_tol=1e-9)
    # at x0 the sigmoid is at half height
    assert math.isclose(pdockq_from_x(DEFAULT.x0), DEFAULT.L / 2 + DEFAULT.b, rel_tol=1e-9)


def test_reference_value_from_paper_regime():
    # pLDDT_if = 90, 100 contacts -> x = 90*log10(100) = 180 -> pDockQ ~ 0.61
    x = 90 * math.log10(100)
    val = pdockq_from_x(x)
    expected = 0.724 / (1 + math.exp(-0.052 * (180 - 152.611))) + 0.018
    assert math.isclose(val, expected, rel_tol=1e-12)
    assert 0.60 < val < 0.63


def test_no_contacts_gives_b():
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[100.0, 0.0, 0.0]])
    out = pdockq_from_arrays(a, np.array([90.0]), b, np.array([90.0]))
    assert out["n_if_contacts"] == 0
    assert out["pdockq"] == DEFAULT.b


def test_synthetic_interface_counts_and_plddt():
    # chain A: 3 residues on a line; chain B: 2 residues, one within 8 A of two A residues,
    # the other far away.  Contacts: (A0,B0) d=5, (A1,B0) d=sqrt(25+9)=5.83 -> 2 contacts.
    a = np.array([[0, 0, 0], [3, 0, 0], [50, 0, 0]], dtype=float)
    b = np.array([[0, 5, 0], [0, 200, 0]], dtype=float)
    pa = np.array([80.0, 60.0, 10.0])
    pb = np.array([70.0, 5.0])
    n, if_plddt = interface_stats(a, pa, b, pb, cutoff=8.0)
    assert n == 2
    # interface residues: A0, A1, B0 -> mean(80, 60, 70) = 70
    assert math.isclose(if_plddt, 70.0)
    out = pdockq_from_arrays(a, pa, b, pb)
    x = 70.0 * math.log10(2)
    assert math.isclose(out["x"], x)
    assert math.isclose(out["pdockq"], pdockq_from_x(x))
    # tiny interface, low x -> close to the floor
    assert out["pdockq"] < 0.03


def test_monotonic_in_contacts_and_plddt():
    v1 = pdockq_from_x(70 * math.log10(10))
    v2 = pdockq_from_x(70 * math.log10(100))
    v3 = pdockq_from_x(90 * math.log10(100))
    assert v1 < v2 < v3
