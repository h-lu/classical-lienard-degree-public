"""Verify the seven return inequalities in the paper using exact integers.

Run with python3 checks/verify_seed.py. Python's standard library suffices.
Taylor approximations are proposals. Full polynomial defects and all decisions
in the verification use exact arithmetic.
"""

from fractions import Fraction as Q
from math import comb, floor, lcm
import sys


SCALE = 2**256
ORDER = 36
LINEAR = tuple(map(Q, (
    "-0.003525444747007079096354237286",
    "-0.003525444747007079096354237285",
)))
HIGHER = tuple(map(Q, (
    "13.8510158", "0.8651558", "-27.5087995", "-1.7814695", "20.3082601",
)))


def ceil_div(a, b):
    return -((-a) // b)


def multiply(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                out[i+j] += a*b
    return out


def taylor_path(coeffs, inv_step, x, y):
    """Construct dyadic coefficients X,Y on the normalized time interval."""
    initial = Q(x, SCALE)
    shifted = []
    for k in range(7):
        ck = sum(coeffs[j]*comb(j, k)*initial**(j-k) for j in range(k, 7))
        shifted.append(floor(SCALE * ck))
    xs, ys = [x], [y]
    powers = [[0] * ORDER for _ in range(7)]
    for n in range(ORDER):
        if n == 0:
            gn = shifted[0]
        else:
            powers[1][n] = xs[n]
            for k in range(2, min(6, n) + 1):
                powers[k][n] = sum(xs[j] * powers[k-1][n-j]
                                   for j in range(1, n-k+2)) // SCALE
            gn = sum(shifted[k] * powers[k][n] for k in range(1, min(6, n)+1)) // SCALE
        divisor = (n+1) * inv_step
        xs.append((ys[n] - gn) // divisor)
        ys.append(-xs[n] // divisor)
    return xs, ys


def residual_norm(numerators, denominator, inv_step, xs, ys):
    """Return N,D with the normalized-time defect norm bounded by N/D."""
    # If p = xs / SCALE, the following Horner polynomial equals
    # denominator * SCALE**6 * G(p), with exact integer coefficients.
    composed = [numerators[6]]
    for j in range(5, -1, -1):
        composed = multiply(xs, composed)
        composed[0] += numerators[j] * SCALE**(6-j)
    factor = denominator * SCALE**5
    ry_norm = 0
    for i in range(ORDER+1):
        dx = (i+1)*xs[i+1] if i < ORDER else 0
        dy = (i+1)*ys[i+1] if i < ORDER else 0
        composed[i] += factor*(inv_step*dx - ys[i])
        ry_norm += factor*abs(inv_step*dy + xs[i])
    return sum(map(abs, composed)) + ry_norm, inv_step*factor*SCALE


def tube(coefficients, error):
    low = coefficients[0] + sum(min(c, 0) for c in coefficients[1:]) - error
    high = coefficients[0] + sum(max(c, 0) for c in coefficients[1:]) + error
    return low, high


def at_time(coefficients, time):
    """Enclose SCALE times a path coordinate at a time in units of 1/SCALE."""
    low = high = coefficients[-1]
    for c in reversed(coefficients[:-1]):
        low = c + (low*time) // SCALE
        high = c + ceil_div(high*time, SCALE)
    return low, high


def event_height(xs, ys, error, max_x, min_speed):
    # The exact path crosses once in this step, with -dx/dt >= min_speed.
    # Find an approximate crossing, then bound its distance to the true event.
    low, high = 0, SCALE
    for _ in range(180):
        mid = (low+high)//2
        xlo, xhi = at_time(xs, mid)
        if xlo > 0:
            low = mid
        elif xhi < 0:
            high = mid
        else:
            break
    x_error = Q(max(abs(xlo), abs(xhi)) + error, SCALE)
    y_error = Q(error, SCALE) + max_x*x_error/min_speed
    ylo, yhi = at_time(ys, mid)
    result = -Q(yhi, SCALE)-y_error, -Q(ylo, SCALE)+y_error
    assert result[0] > 0 and result[1]-result[0] < Q("1e-40")
    return result


def check_derivative(coeffs):
    # G' = c1 + x*q(x). Complete the square twice for q.
    q0, q1, q2, q3, q4 = [j*coeffs[j] for j in range(2, 7)]
    u, v = q3/(2*q4), Q(-23, 50)
    w = q2 - q4*(u*u+2*v)
    z = (q1-2*q4*u*v)/(2*w)
    remainder = q0-q4*v*v-w*z*z
    assert q4 > 0 and w > 0 and remainder > 1
    q = [q0, q1, q2, q3, q4]
    assert abs(coeffs[1]) + sum(abs(c)/100**(j+1) for j, c in enumerate(q)) < 1


def return_height(coeffs, x0, height):
    check_derivative(coeffs)
    denominator = lcm(*(c.denominator for c in coeffs))
    numerators = [int(c * denominator) for c in coeffs]
    x, y = floor(x0*SCALE), floor(height*SCALE)
    error = 2  # sum of both initial rounding errors in units of 1/SCALE
    elapsed, steps, inv_step = Q(0), 0, 16
    while elapsed < 40:
        assert inv_step < 2**32
        xs, ys = taylor_path(coeffs, inv_step, x, y)
        num, den = residual_norm(numerators, denominator, inv_step, xs, ys)
        if num * inv_step * 10**53 >= den:
            inv_step *= 2
            continue
        # e^h <= 1/(1-h). The vector field's logarithmic norm is <= 1.
        local_error = ceil_div(SCALE*num, den)
        end_error = ceil_div((error+local_error)*inv_step, inv_step-1)
        xmin, xmax = tube(xs, end_error)
        if 100*xmin <= -SCALE:
            inv_step *= 2
            continue
        xend, yend = sum(xs), sum(ys)
        if xmin <= 0:
            ymin, ymax = tube(ys, end_error)
            max_x = Q(max(abs(xmin), abs(xmax)), SCALE)
            max_g = sum(abs(c)*max_x**j for j, c in enumerate(coeffs))
            departing = steps == 0 and x0 == 0 and Q(ymin, SCALE) > max_g
            min_speed = -Q(ymax, SCALE) - max_g
            if not departing:
                if min_speed <= 0 or x <= error:
                    inv_step *= 2
                    continue
                if xend < -end_error:
                    result = event_height(xs, ys, end_error, max_x, min_speed)
                    return result, steps+1
                if xend <= end_error:
                    inv_step *= 2
                    continue
        x, y, error = xend, yend, end_error
        elapsed += Q(1, inv_step)
        steps += 1
        if inv_step > 8 and num * inv_step * 10**66 < den:
            inv_step //= 2
    raise ArithmeticError("return time exceeded 40")


def verify(label, linear, height, plus_x, minus_x, lower, upper):
    coeffs = (Q(0), linear, *HIGHER)
    mirror = tuple((-1)**j*c for j, c in enumerate(coeffs))
    plus, psteps = return_height(coeffs, Q(plus_x), Q(height))
    minus, msteps = return_height(mirror, Q(minus_x), Q(height))
    lo, hi = plus[0]-minus[1], plus[1]-minus[0]
    assert Q(lower) < lo <= hi < Q(upper), label
    print(f"{label:<21} ({lower}, {upper})   {psteps} + {msteps} steps", flush=True)


def main():
    low, high = LINEAR
    assert high-low == Q("1e-30")
    # Each row bounds A_plus(plus_x,height) - A_minus(minus_x,height).
    checks = [
        ("h0, upper parameter", high, "0.5", 0, 0, "2e-7", "3e-7"),
        ("h1, lower parameter", low, "3.7", 0, 0, "-3e-24", "-1e-24"),
        ("h2, upper parameter", high, "3.88", 0, 0, "3e-27", "4e-27"),
        ("h3, lower parameter", low, "4.03", 0, 0, "-3e-28", "-2e-28"),
        ("h4, upper parameter", high, "4.36", 0, 0, "6e-32", "8e-32"),
        ("high return, lower", low, "8", 0, 2, "2e-32", "3e-32"),
        ("high return, upper", high, "8", 2, 0, "-8e-32", "-6e-32"),
    ]
    if not __debug__:
        sys.exit("Run Python with assertions enabled.")
    for row in checks:
        verify(*row)


if __name__ == "__main__":
    main()
