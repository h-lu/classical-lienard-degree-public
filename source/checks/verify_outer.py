"""Check the moving-endpoint action at three octic samples using rational bounds.

The quantity checked is the spatial integral on [R^-2,Y] in the proof.
The trajectory-to-integral estimate is proved in Proposition 3 of the paper.
"""

from fractions import Fraction as Q
import sys

from verify_seed import HIGHER, LINEAR


R = 10_000
LEFT, RIGHT = Q(1, 100), Q(2)
A2, A3, A4, A5, A6 = HIGHER
E8 = Q(15, 2)*A6


def value(coeffs, z):
    result = Q(0)
    for c in reversed(coeffs):
        result = result*z+c
    return result


def inverse_bracket(coeffs, height):
    lo, hi = LEFT, RIGHT
    assert value(coeffs, lo) < height < value(coeffs, hi)
    for _ in range(110):
        mid = (lo+hi)/2
        if value(coeffs, mid) < height:
            lo = mid
        else:
            hi = mid
    return lo, hi


def action(sign, height):
    a = LINEAR[0]
    phi = [Q(0), sign*a/R**5, A2/R**4, sign*A3/R**3, A4/R**2,
           sign*A5/R, A6, sign*A5/R, E8]
    # A positive lower bound for Phi'(z)/z^5 on [LEFT,RIGHT].
    rz = R*LEFT
    lower_slope = (6*A6 - abs(A5)*(5/LEFT+7*RIGHT)/R
                   - 4*abs(A4)/rz**2 - 3*abs(A3)/rz**3
                   - 2*abs(A2)/rz**4 - abs(a)/rz**5)
    assert lower_slope > 0
    z0lo, z0hi = inverse_bracket(phi, Q(1,R**2))
    z1lo, z1hi = inverse_bracket(phi, height)
    assert z0hi < z1lo
    derivative = [k*c for k, c in enumerate(phi) if k]
    square = [Q(0)]*(2*len(derivative)-1)
    for i, b in enumerate(derivative):
        for j, c in enumerate(derivative):
            square[i+j] += b*c
    polynomial = [Q(0), *(c/k for k,c in enumerate(square) if k)]
    # The exact primitive is polynomial(z) + square[0]*log(z).
    # Its derivative is nonnegative. For both rational endpoint pairs,
    # 0 <= log(z1/z0) <= RIGHT/LEFT-1.
    low = value(polynomial,z1lo)-value(polynomial,z0hi)
    high = value(polynomial,z1hi)-value(polynomial,z0lo)+square[0]*(RIGHT/LEFT-1)
    return low, high


def main():
    if not __debug__:
        sys.exit("Run Python with assertions enabled.")
    for x, expected_sign in [(Q(1,4),-1), (Q(3,5),1), (Q(1),-1)]:
        height = A6*x**6+E8*x**8
        w = 6*A6*A5*x**9*(Q(1,9)-7*x**2/11+10*x**4/13)
        assert expected_sign*w > 0
        plo, phi = action(1,height)
        mlo, mhi = action(-1,height)
        normalized = sorted((R*(plo-mhi)/(2*w), R*(phi-mlo)/(2*w)))
        assert Q(9999,10_000) < normalized[0] <= normalized[1] < Q(10001,10_000)
        print(f"X={x}, sign {expected_sign:+d}, normalized action in (0.9999, 1.0001)")


if __name__ == "__main__":
    main()
