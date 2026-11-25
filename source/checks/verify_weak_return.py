"""Exact rational checks for the weak full-return obstruction and its degree-12 witness."""
from fractions import Fraction as Q
from math import comb, factorial
import sys

from verify_seed import multiply


def mean_sine_power(k):
    """Mean of sin(theta)^(2k+2), with the full period factored out."""
    return Q(comb(2*k+2, k+1), 4**(k+1))


def coefficient_radius(n):
    d = (n-1)//2
    rho = Q(1, 2**(n+6))
    determinant_factor = factorial(d+1)*(d+1)*2**d
    gap = Q(1, n**(d*(d+1)//2))
    threshold = rho**2*gap/(72*(n+2)*factorial(d)*determinant_factor)
    return rho, threshold


def check_wronskian_bounds(n):
    d = (n-1)//2
    rho, threshold = coefficient_radius(n)
    entry_error = 36*(n+2)*factorial(d)*threshold/rho**2
    assert 3*threshold < rho/4 and entry_error < 1
    sine_mean = Q(1)
    determinant = Q(1)
    for k in range(d+1):
        sine_mean *= Q(2*k+1, 2*k+2)
        assert sine_mean == mean_sine_power(k)
        # Direct rational Vandermonde product for the first k+1 columns.
        for i in range(k):
            determinant *= Q(2*(k-i), n)
        error = factorial(k+1)*(k+1)*(1+entry_error)**k*entry_error
        assert error <= determinant/2
        assert 6*sine_mean > Q(1, n+2)


def evaluate(coeffs, t):
    result = Q(0)
    for c in reversed(coeffs):
        result = result*t+c
    return result


def finite_return_error(n, norm, odd_norm_scaled):
    """Error for (P(h)-h)/(h*tau^5) about the signed sine moment.

    Differentiate the true radial ODE along the segment from the even
    polynomial to the full polynomial. The odd coefficient norm here has
    already been divided by tau^5. See the verification record.
    """
    qbound = 2**(n-1)*norm
    assert qbound <= Q(1, 128)
    denominator = 1-qbound
    radius_drift = 16*qbound/denominator
    radial_derivative = n*qbound/denominator**2
    assert radius_drift < Q(1, 2) and 8*radial_derivative < 1
    weight_error = 8*radial_derivative/(1-8*radial_derivative)
    odd_bound = 2**(n-1)*odd_norm_scaled
    integrand_error = n*odd_bound*radius_drift+odd_norm_scaled*(2*qbound+qbound**2)
    return 8*(integrand_error+2*odd_bound*weight_error)/denominator**2


def main():
    if not __debug__:
        sys.exit("Run with assertions enabled.")
    for n in (1, 2, 7, 12, 16):
        check_wronskian_bounds(n)

    n, tau = 12, Q(1, 10**60)
    d = (n-1)//2
    centers = [Q(i, 6) for i in range(1, d+1)]
    target = [Q(1)]
    for h in centers:
        target = multiply(target, [-h*h, Q(1)])
    odd = [a/mean_sine_power(k) for k, a in enumerate(target)]
    coefficients = [Q(0)]*(n+1)
    coefficients[n] = tau
    for k, c in enumerate(odd):
        coefficients[2*k+1] = tau**5*c
    rho, threshold = coefficient_radius(n)
    full_norm = sum(abs(c) for c in coefficients)
    assert coefficients[n] and coefficients[n-1] and full_norm < threshold
    assert all(evaluate(target, h*h) == 0 for h in centers)
    interpolation_bound = Q((4*n)**n, factorial(n))
    damping_threshold = threshold/(2**(n-1)*interpolation_bound)
    box_damping = sum(j*abs(coefficients[j])*2**(j-1) for j in range(1, n+1))
    assert box_damping < damping_threshold

    signs = []
    margins = []
    for h0 in centers:
        row = []
        for h in (h0-Q(1, 24), h0+Q(1, 24)):
            odd_norm = sum(abs(c)*h**(2*k) for k, c in enumerate(odd))
            norm = sum(abs(c)*h**(j-1) for j, c in enumerate(coefficients) if c)
            assert norm == tau*h**11+tau**5*odd_norm
            q = evaluate(target, h*h)
            # Two separate enclosures, from complex bounds and the real variation equation.
            complex_error = 4*norm*odd_norm/rho**2
            real_error = finite_return_error(n, norm, odd_norm)
            assert complex_error < Q(4, 10**49) and real_error < Q(2, 10**50)
            # -2*pi*q lies between -6*q and -8*q in either order.
            center_lo, center_hi = sorted((-6*q, -8*q))
            endpoint_signs = []
            for error in (complex_error, real_error):
                lo, hi = center_lo-error, center_hi+error
                assert hi < 0 or lo > 0
                margins.append(min(abs(lo), abs(hi)))
                endpoint_signs.append(1 if lo > 0 else -1)
            assert endpoint_signs[0] == endpoint_signs[1]
            row.append(endpoint_signs[0])
        signs.append(row)
    assert signs == [[1, -1], [-1, 1], [1, -1], [-1, 1], [1, -1]]
    assert min(margins) > Q(3, 10**4)
    assert d == 5
    print("Sine moments and all prefix determinant margins pass at n = 1, 2, 7, 12, 16.")
    print("The rational degree-12 polynomial is inside the explicit weak radius for all 0 < h <= 1.")
    print("Physical damping on |x| <= 2 is below the intrinsic orbit threshold.")
    print("Complex and real full-return bounds certify five disjoint brackets, each with one simple zero.")
    print("The specified weak-damping polynomial has five hyperbolic cycles in the stated height range.")


if __name__ == '__main__':
    main()
