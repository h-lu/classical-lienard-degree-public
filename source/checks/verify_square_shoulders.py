"""Exact checks for the square-derivative moving-shoulder construction."""

from fractions import Fraction as Q
from math import ceil, comb, floor
import sys

from verify_outer import value
from verify_seed import HIGHER, LINEAR, multiply


BASE = 25
SAMPLES = tuple(map(Q, ("-12", "-7.62", "-7.61", "0", "9.99", "10.01", "100", "100000")))
WIDTHS = ((81, 82), ("-.85", "-.82"), (".40", ".44"), (-79, -77),
          (21, 23), ("-.43", "-.40"), (610, 615), (-98000, -96000))


def pair(a, b):
    return (a, b) if a <= b else (b, a)


def plus(a, b):
    return a[0]+b[0], a[1]+b[1]


def minus(a, b):
    return a[0]-b[1], a[1]-b[0]


def times(a, b):
    ends = [x*y for x in a for y in b]
    return min(ends), max(ends)


def inside(interval, lo, hi):
    assert Q(lo) < interval[0] <= interval[1] < Q(hi), (interval, lo, hi)


def root(f, a, b, steps=110):
    """Bisect a strict sign bracket. The function may decrease."""
    fa, fb = f(a), f(b)
    assert fa*fb < 0
    for _ in range(steps):
        c = (a+b)/2
        fc = f(c)
        if not fc:
            return c, c
        if fc*fa > 0:
            a = c
        else:
            b = c
    return a, b


def cubic_rows(a, shift, height):
    k = root(lambda x: 3*x*x+a, Q(0), 1+abs(a))
    min_height = plus((shift, shift), times((-2, -2), times(k, times(k, k))))
    max_height = minus((2*shift, 2*shift), min_height)
    assert all(height < v[0] or height > v[1] for v in (min_height, max_height))
    mid = min_height[1] < height < max_height[0]
    far = 2+abs(a)
    while far**3+a*far <= abs(height-shift):
        far *= 2
    if height < min_height[0]:
        bounds = (-far, -k[1])
    else:
        bounds = (k[1], far)
    v = root(lambda x: x**3+a*x+shift-height, *bounds)

    def primitive(x):
        return Q(9,5)*x**5+2*a*x**3+a*a*x

    def at(x):
        # The derivative of the primitive is a square.
        return primitive(x[0]), primitive(x[1])

    landing = at(times((-2, -2), k))
    incoming = at(v)
    if height > min_height[1]:
        incoming = minus(incoming, minus(at(k), landing))
    return (landing if mid else incoming), incoming


def cubic_certificate():
    up, inc = [], []
    for height, band in zip(SAMPLES, WIDTHS):
        p, n = cubic_rows(-1, Q(-8), height), cubic_rows(-3, Q(8), height)
        u, v = [minus(minus(x,y), (Q(493,2),)*2) for x,y in zip(p,n)]
        inside(u, *band)
        if height in (SAMPLES[i] for i in (0, 3, 6, 7)):
            assert v == u
        print(f"split cubic xi={height}: upper action in ({band[0]}, {band[1]})")
        assert v[0]*v[1] > 0
        up.append(u[0] > 0)
        inc.append(v[0] > 0)
    assert up == [True, False]*4
    assert [inc[i] for i in (0,3,6,7)] == [True, False, True, False]
    # A second, much larger skew exercises the moving choice of mean slope.
    c, alpha, level = Q(19683), Q("-9.588"), Q(82385566)
    samples = list(map(Q, ("-100000","-19673.32","-19673.31","0","19696.25",
                           "19696.27",str(c**2),"10000000000000000")))
    for i, t in enumerate(samples):
        p, n = cubic_rows(alpha+1,-c,t), cubic_rows(alpha-1,c,t)
        for k in range(2):
            if k == 0 or i in (0,3,6,7):
                bounds = minus(minus(p[k],n[k]), (level, level))
                assert min((-1)**i*v for v in bounds) > 0
    print("Rational cubic samples at skews 8 and 19683 give seven upper and three incoming changes.")


def linear_combination(weights, polys):
    out = [Q(0)]*max(map(len, polys))
    for w, poly in zip(weights, polys):
        for i, a in enumerate(poly):
            out[i] += w*a
    return out


def monomial(n):
    return [Q(0)]*n+[Q(1)]


def deriv(poly):
    return [i*v for i,v in enumerate(poly) if i]


def primitive(poly):
    return [Q(0)]+[v/Q(i+1) for i,v in enumerate(poly)]


def moment(f, g, z):
    """Integral of f dg from zero to z."""
    return value(primitive(multiply(f, deriv(g))), z)


def solve(rows, rhs):
    data = [list(row)+[b] for row,b in zip(rows,rhs)]
    n = len(data)
    for i in range(n):
        k = next(k for k in range(i, n) if data[k][i])
        data[i], data[k] = data[k], data[i]
        pivot = data[i][i]
        data[i] = [v/pivot for v in data[i]]
        for k in range(n):
            if k != i:
                weight = data[k][i]
                data[k] = [v-weight*w for v,w in zip(data[k],data[i])]
    return [row[-1] for row in data]


def profiles(shoulders):
    b = [v*Q(3,2) for v in monomial(BASE)]
    for q in shoulders:
        b = multiply(b, [1, 0, -1/q**2])
    a = multiply(b,b)
    even = primitive([Q(0)]+a)
    odd_modes = [monomial(BASE+2+2*i) for i in range(len(shoulders)+1)]
    rows = [[moment(a, mode, q) for mode in odd_modes] for q in shoulders]
    coeff = [Q(-1), *solve([r[1:] for r in rows], [r[0] for r in rows])]
    h = linear_combination(coeff, odd_modes)
    assert all(moment(a,h,q) == 0 for q in shoulders)
    assert all(value(h,q)*value(deriv(b),q) < 0 for q in shoulders)
    assert h[-1]*b[-1] < 0
    return b, a, even, h


def control(b, shoulders, values, moments):
    # The lowest z^(2j0) coefficient is fixed at zero to remove the square kernel.
    modes = [monomial(2*BASE+2*i) for i in range(1, 2*len(shoulders)+1)]
    rows = [[value(p,q) for p in modes] for q in shoulders]
    rows += [[moment(p,b,q) for p in modes] for q in shoulders]
    out = linear_combination(solve(rows, [*values,*moments]), modes)
    assert [value(out,q) for q in shoulders] == values
    assert [moment(out,b,q) for q in shoulders] == moments
    return out


def flight_integrals(shoulders, even, b, radii, expected_signs):
    # Normalize by w^2/r^(BASE+2) before adding panel bounds on a dyadic lattice.
    j = len(b)-1
    zmax = 8*shoulders[-1]
    pole_product = Q(1)
    for q in shoulders:
        pole_product *= q*q
    leading = Q(3,2)/pole_product
    lower = leading**2*(1-Q(1,2**(2*j+2)))/((2*j+2)*4**len(shoulders))
    tail = leading*zmax**(-3*j-2)/(lower**2*(3*j+2))
    grid = 2**100
    ends = [Q(0)]
    z = Q(1,128)
    ends.append(z)
    while z < zmax:
        ends.extend(z*(1+Q(i,32)) for i in range(1,33))
        z *= 2
    assert ends[-1] == zmax
    energies = [value(even,z) for z in ends]
    margins = []
    for r, expected in zip(radii, expected_signs):
        w = value(even,r)
        norm = w*w/r**(BASE+2)
        low, high = 0, 0
        for i, (a,c) in enumerate(zip(ends,ends[1:])):
            num = Q(3,2)*a**(BASE+1), Q(3,2)*c**(BASE+1)
            for q in shoulders:
                num = times(num,(1-(c/q)**2,1-(a/q)**2))
            part = times(num, (norm*(c-a)/(w+energies[i+1])**2,
                               norm*(c-a)/(w+energies[i])**2))
            low += floor(grid*part[0])
            high += ceil(grid*part[1])
        margin = min(expected*Q(v, grid) for v in (low, high))/norm - tail
        assert margin > 0
        margins.append((w, margin))
    return margins


def physical_flights(b, even, h, gamma, d, zcontrol):
    """Bound the full negative-flight Volterra remainder with the sextic seed."""
    radii = (Q(1,2), Q(3))
    margins = flight_integrals([Q(1)], even, b, radii, (1,-1))
    eta = Q(1, 10**600)
    sqrt_eta = Q(1, 10**300)
    seed_scale = Q(1, 10**960)  # eta^(8/5)
    space_scale = eta**-10
    eps = seed_scale**2/space_scale**10
    m = len(even)-1
    assert eta == sqrt_eta**2 and seed_scale**5 == eta**8 and m == 56

    # For z>=4 integrate on [z/2,z], where |B(v)| >= (3/4)v^27.
    tail_constant = Q(9,16)*(1-Q(1, 2**m))/m
    coefficient_norm = lambda poly: sum(map(abs, poly))
    even_error = (eta*coefficient_norm(primitive([Q(0)]+d))
                  + eta*sqrt_eta*coefficient_norm(primitive([Q(0)]+zcontrol))
                  + seed_scale*coefficient_norm(HIGHER[::2]))
    odd_error = (eta*sqrt_eta*gamma*coefficient_norm(h)
                 + seed_scale/space_scale*(coefficient_norm(HIGHER[1::2])
                                           + max(map(abs, LINEAR))))
    for w, margin in margins:
        lower = min(w/(1+4**m), tail_constant/2)
        # w+E >= lower*(1+z^m), and integral z/(1+z^m) < 1.
        raw_relative = (even_error + eta*coefficient_norm(b) + odd_error)/lower
        kinetic_relative = 2*eps/lower**2
        relative = raw_relative + kinetic_relative
        assert relative < Q(1,2) and 4*eps/lower**2 < 1
        # Expand (1+x)^(-1)=1-x+x^2/(1+x) on both sides.
        # The common even part cancels and |x^2/(1+x)| <= 2*relative^2.
        error = (odd_error/lower + kinetic_relative + 2*relative**2)/(lower*eta)
        assert error < margin
    print("Full seeded half-flights at eta=10^-600 have terminal signs -,+ uniformly on the shooting interval.")


def poly_box(poly, interval):
    result = (Q(0), Q(0))
    for coefficient in reversed(poly):
        result = plus((coefficient, coefficient), times(result, interval))
    return result


def finite_actions(b, a, h):
    # At z=1 the u-derivative of B is -3, with l_i=1, d_i=3 and r_i=9.
    q = Q(1)
    gamma = 24/value(h,q)
    assert gamma > 0 and value(deriv(b),q) == -3
    d = control(b, [q], [Q(-6)], [Q(0)])
    target = Q(9*493,8)-gamma*moment(d,h,q)
    zcontrol = control(b, [q], [Q(0)], [target])
    # No sextic seed in this finite singular-action test.
    lam = Q(1,10**12)
    eta = lam**2
    slope = linear_combination([1,eta,eta*lam], [a,d,zcontrol])
    even = primitive([Q(0)]+slope)
    odd = linear_combination([1,gamma*lam], [b,h])
    center = value(even,q)
    tol = Q(1,10**8)

    def branch_data(s):
        profile = linear_combination([1,s*eta], [even,odd])
        derivative = deriv(profile)
        primitive_action = primitive(multiply(derivative,derivative)[1:])
        # In v=(z-1)/sqrt(eta) the three root intervals stay bounded.
        deriv_v = shift_at_one(derivative,lam)
        prof_v = shift_at_one(profile,lam)
        act_v = shift_at_one(primitive_action,lam)
        a_cubic = -2-s  # physical s=+ has negative normalized shoulder slope
        discrepancy = linear_combination([1/(3*eta), -1], [deriv_v, [a_cubic,0,3]])
        for polynomial in (discrepancy, deriv(discrepancy)):
            assert sum(abs(c)*70**i for i,c in enumerate(polynomial)) < Q(1,100)
        # The normalized slope is negative on [-.3,.3] and strictly
        # monotone on each of [-70,-.3] and [.3,70].
        k0 = Q(1) if a_cubic == -3 else Q(".57735")
        loq, hiq = root(lambda v:value(deriv_v,v), k0-Q(".02"), k0+Q(".02"), 140)
        mu = poly_box(prof_v,(loq,hiq))
        min_t = pair((mu[0]-center)/(3*eta*lam), (mu[1]-center)/(3*eta*lam))
        max_v = root(lambda v:value(deriv_v,v), -k0-Q(".02"), -k0+Q(".02"),140)
        peak = poly_box(prof_v,max_v)
        max_t = pair((peak[0]-center)/(3*eta*lam),(peak[1]-center)/(3*eta*lam))
        r1 = root(lambda v:value(prof_v,v)-mu[0], -2*k0-Q(".03"), -2*k0+Q(".03"),140)
        r2 = root(lambda v:value(prof_v,v)-mu[1], -2*k0-Q(".03"), -2*k0+Q(".03"),140)
        land = poly_box(act_v,(r1[0],r2[1]))
        skip = minus(poly_box(act_v,(loq,hiq)),land)

        def actions_at(height):
            assert all(height < v[0] or height > v[1] for v in (min_t, max_t))
            if height < min_t[0]:
                start, end = Q(-70), max_v[0]
            else:
                start, end = hiq, Q(70)
            actual = root(lambda v:value(prof_v,v)-center-3*eta*lam*height, start,end,140)
            incoming = poly_box(act_v, actual)
            if height > min_t[1]:
                incoming = minus(incoming, skip)
            return (land if min_t[1] < height < max_t[0] else incoming), incoming

        return actions_at

    positive, negative = branch_data(1), branch_data(-1)
    for i,t in enumerate(SAMPLES):
        p, n = positive(t), negative(t)
        for k in range(2):
            d_scaled = times((1/(9*lam**5),)*2, minus(n[k],p[k]))
            if k==0:
                inside(d_scaled,*WIDTHS[i])
            elif i in (0,3,6,7):
                assert min((1 if i in (0,6) else -1)*v for v in d_scaled) > tol
    print("Finite polynomial singular actions give the seven and three changes.")
    return gamma, d, zcontrol


def shift_at_one(poly, scale):
    return [sum(v*comb(i,k) for i,v in enumerate(poly) if i>=k)*scale**k
            for k in range(len(poly))]


def degrees():
    for h in (1,2,12,119,200):
        m, j = 52+4*h, 27+2*h
        k = (m-j-1)//2
        q, r = 7*h+6, 4*h+3
        upper = [1,-1,1,-1,1]
        incoming = []
        for i in range(1, h+1):
            orientation = (-1)**i
            upper.extend(orientation*(-1)**t for t in range(8))
            incoming.extend(-orientation*(-1)**t for t in range(4))
        upper.append(-upper[-1])
        incoming.reverse()
        incoming.extend([-1,1])
        incoming.extend(-(-1)**i for i in range(h+1))
        assert sum(x != y for x,y in zip(upper,upper[1:])) == q
        assert sum(x != y for x,y in zip(incoming,incoming[1:])) == r
        assert incoming[-1] == -(-1)**h  # opposite to the leading sign of B
        assert m-j == 2*k+1 and k == h+12
        assert q+r+k+1 == 12*h+22
        assert q+r+k+1+k+2+(k%2) == 13*h+36+(h%2)
        assert m+5 == 57+4*h
    assert 13*119+36+1 == 1584 and 12*119+154 == 1582
    for e, offset in enumerate((154,157,158,161)):
        first = next(h for h in range(1,130) if 13*h+36+h%2 > 12*h+offset)
        assert first == (119,121,123,125)[e]
        # Raising h by two preserves parity and increases the surplus by two.
        for h in (125,126):
            assert 13*h+36+h%2 > 12*h+offset
            assert 4*(13*h+36+h%2) >= 13*(57+4*h+e)-636
    print("Degree and both transfer counts agree, including parity and terminal sign.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run with assertions enabled.")
    cubic_certificate()
    shoulders = [Q(1),Q(16)]
    b, a, e, h = profiles(shoulders)
    control(b, shoulders, [Q(-6), Q(2)], [Q(1), Q(-3)])
    flight_integrals(shoulders, e, b, (Q(1,2), Q(3), Q(32)), (1,-1,1))
    print("Two-shoulder flight integrals have the whole-ray signs +,-,+.")
    b, a, e, h = profiles([Q(1)])
    physical_flights(b, e, h, *finite_actions(b,a,h))
    degrees()
