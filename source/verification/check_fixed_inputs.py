"""Reconstruct the fixed fields and verify escape at the rational enclosure endpoint."""

from fractions import Fraction as Q
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "checks"))
from verify_joint_shoulder_cap import finite_profile, HIGHER, LINEAR, linear_combination
from verify_outer import value
from verify_remote_odd_powers import bernstein_positive, shifted


def check_outer_enclosure(data, scaled, d, lam, height, epsilon):
    """Check the forward comparison in section 5 of fixed-degree12/PROOF.md."""
    eta, rho = lam**2, lam**4
    at_d = value(scaled, d)
    cap_bound = at_d + 5*10**10*eta**2
    lo, hi = (Q(h)/height for h in data["outer_upper_axis_interval"])
    assert (lo, hi) == (at_d-rho, cap_bound+2*rho)
    assert at_d > Q("0.012")

    # Positive Bernstein coefficients prove P_+(z) < cap_bound on [0,2].
    cap_margin = [-c for c in scaled]
    cap_margin[0] += cap_bound
    grid = (Q(0), Q(1, 2), Q(1), d-100*eta, d+100*eta, Q(2))
    assert all(bernstein_positive(cap_margin, a, b) for a, b in zip(grid, grid[1:]))
    assert 2*epsilon/rho < rho

    # All coefficients in powers of z-2 are negative after adding the tail bound.
    tail = scaled.copy()
    tail[-1] += Q(1, 1000)
    assert all(c < 0 for c in shifted(tail, Q(2)))
    tail_loss = epsilon / (Q(1, 1000)*2**6*5*2**5)
    assert tail_loss == Q(25, 256)*epsilon < rho
    print("Rational enclosure: both endpoints agree with the formulas in scaled height.")
    print("Plus cap: five Bernstein intervals certify M_+ < P_+(d) + m_lift on 0 <= z <= 2.")
    print("Plus escape at the rational upper endpoint: height losses <= 2*epsilon/rho < rho "
          "and <= 25*epsilon/256 < rho. The half-return is undefined there.")


def main():
    for name, orientation in (("outer-variant", -1), ("fixed-diagnostic", 1)):
        path = HERE / "fixed-degree12" / "inputs" / f"{name}-field.json"
        data = json.loads(path.read_text())
        d, theta, lam = (Q(data[key]) for key in ("d", "theta", "lambda"))
        assert (d, theta, lam) == (Q(6, 5), orientation*Q("1e-20"), Q("1e-40"))
        even, odd, _, cutoff = finite_profile(d, theta, lam)
        length = 20 / cutoff
        height = 6*HIGHER[-1]*length**6
        scaled = linear_combination([1, 1], [even, odd])
        physical = [height*c / length**i for i, c in enumerate(scaled)]
        epsilon = length**2/height**2
        assert (length, height, epsilon) == tuple(Q(data[k]) for k in ("L", "S", "epsilon"))
        assert physical == list(map(Q, data["coefficients"]))
        assert len(physical) == 13 and physical[-1] < 0
        assert physical[:2] == [0, sum(LINEAR)/2]
        print(f"{name}: all 13 coefficients match the Section S5 generator; theta={orientation:+d}e-20, degree 12.")
        if orientation < 0:
            check_outer_enclosure(data, scaled, d, lam, height, epsilon)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()
