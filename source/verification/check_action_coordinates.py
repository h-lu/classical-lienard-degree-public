"""Exact algebra checks for the actions in Proposition 3.3 and Sections S3 and S5."""

from fractions import Fraction as Q
import json
from pathlib import Path
import sys

FIXED = Path(__file__).resolve().parent / "fixed-degree12"
sys.path.insert(0, str(FIXED))
from rational import at, deriv, mul, primitive, sub


def integral(poly, left, right):
    antiderivative = primitive(poly)
    return at(antiderivative, right) - at(antiderivative, left)


def even_action(filename):
    data = json.loads((FIXED / "inputs" / filename).read_text())
    length, height = Q(data["L"]), Q(data["S"])
    phi = [Q(c) * length**k / height if k % 2 == 0 else Q(0)
           for k, c in enumerate(data["coefficients"])]
    p = [c * 2**k for k, c in enumerate(phi[::2])]
    dz, du = deriv(phi), deriv(p)
    square = mul(dz, dz)
    assert square[0] == 0
    spatial_primitive = primitive(square[1:])  # density Phi'(z)^2/z
    u_primitive = primitive(mul(du, du))
    # Compare the polynomials after substitution, with both constants fixed at zero.
    assert all(c == 0 for c in spatial_primitive[1::2])
    assert spatial_primitive[::2] == [Q(c) / 2**k for k, c in enumerate(u_primitive)]

    left, right = Q(6, 5), Q(3, 2)  # test endpoints in z
    spatial = at(spatial_primitive, right) - at(spatial_primitive, left)
    assert spatial == at(u_primitive, right**2 / 2) - at(u_primitive, left**2 / 2)
    assert spatial != at(u_primitive, right) - at(u_primitive, left)
    for z, sign in ((Q(1, 4), 1), (right, -1)):
        slope = at(dz, z) / z
        assert slope == at(du, z**2 / 2) and sign * slope > 0


def first_variation():
    # Limiting degree-twelve profile in (S5.1), with a test odd direction z^5.
    d, v = Q(6, 5), Q(3, 2)
    a = [0] * 4 + mul([1, 0, -2, 0, 1], [1, 0, -1 / d**2])
    b = [0] * 5 + [1]
    ab_prime, a_prime_b = mul(a, deriv(b)), mul(deriv(a), b)
    endpoint_velocity = -at(b, v) / (v * at(a, v))
    endpoint_term = endpoint_velocity * v * at(a, v)**2
    coefficient = 2 * integral(ab_prime, d, v) + endpoint_term
    assert at(a, d) == 0 and endpoint_term != 0
    assert coefficient == integral(sub(ab_prime, a_prime_b), d, v)
    assert coefficient == -integral(sub(a_prime_b, ab_prime), d, v) != 0

    # Direct differentiation of DMH's right-to-left integral for F=z^6+delta*z^5.
    # At fixed right endpoint 1, the equal-height left endpoint moves at -2/6.
    moving_left = Q(-2, 6)
    signed_derivative = integral([0] * 8 + [60], 1, -1) - 36 * moving_left
    w = integral([0] * 8 + [6], 0, 1)
    assert signed_derivative == -2 * w == Q(-4, 3)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    for filename in ("fixed-diagnostic-field.json", "outer-variant-field.json"):
        even_action(filename)
        print(f"{filename}: the two action primitives agree in z. Mixed endpoints give a different value.")
    first_variation()
    print("Falling-interval endpoint variation equals minus the DMH primitive increment.")
    print("For F=z^6+delta*z^5 at x=1, dI/ddelta at delta=0 is -4/3 and W(1)=2/3.")
