"""Exact spatial and degree checks for a moving cap coupled to a cubic shoulder."""

from fractions import Fraction as Q
import sys

from verify_entrance import SHIFT, TAIL
from verify_outer import value
from verify_remote_odd_powers import bernstein_positive, polynomial
from verify_seed import HIGHER, LINEAR, multiply, return_height
from verify_square_shoulders import (
    SAMPLES, WIDTHS, deriv, inside, linear_combination, minus, moment, monomial,
    plus, poly_box, primitive, root, shift_at_one, times,
)
from verify_unequal_cubic_layers import OFFSET, outward_signs, scaled, sign


ONE = Q(1)
SHOULDER = polynomial({0: 1, 2: -1})
FACES = (Q(23, 20), Q(5, 4))
TIGHT = (Q('1.19888999543'), Q('1.19888999544'))


def difference(left, right):
    return linear_combination([1, -1], [left, right])


def integral(p, x):
    return value(primitive(p), x)


def wronskian(a, b):
    return difference(multiply(a, deriv(b)), multiply(deriv(a), b))


def profile(d, theta=Q(0)):
    cap = polynomial({0: 1, 2: -1/d**2})
    both = multiply(SHOULDER, cap)
    a = [0]*4 + multiply(SHOULDER, both)
    e = primitive([0]+a)
    omega = [0]*8 + multiply(both, both)
    base = multiply(omega, polynomial({0: 1, 2: 1}))
    tilt = multiply(omega, polynomial({2: 3, 4: -1}))
    k = -integral(base, 1)/integral(tilt, 1)
    w = linear_combination([1, k], [base, tilt])
    seed = scaled(-1, [0]*5 + multiply(both, polynomial({0: 1, 2: k})))
    assert k < -1 and wronskian(a, seed) == scaled(-1, w)
    u = [0]*7 + SHOULDER
    u1 = multiply(u, cap)
    assert moment(a, u1, 1) > 0
    n = scaled(1/value(u, d), linear_combination(
        [1, -moment(a, u, 1)/moment(a, u1, 1)], [u, u1]))
    assert value(n, 1) == moment(a, n, 1) == 0 and value(n, d) == 1
    unscaled = linear_combination([1, theta], [seed, n])
    unit = 4*(1-1/d**2)/3
    beta = -unit/value(deriv(unscaled), 1)
    b = scaled(beta, unscaled)
    assert beta > 0 and b[5] == -beta and value(b, d) == beta*theta
    assert value(b, 1) == moment(a, b, 1) == 0 and value(deriv(b), 1) == -unit
    h0 = scaled(8*unit, monomial(7))
    h = linear_combination([1, -moment(a, h0, 1)/moment(a, u1, 1)], [h0, u1])
    even_mode = b[1:]
    pivot = moment(even_mode, b, 1)
    assert pivot == integral(multiply(even_mode, even_mode), 1)/2 > 0
    d0 = scaled(-2*unit, monomial(4))
    control = linear_combination([1, -moment(d0, b, 1)/pivot], [d0, even_mode])
    z_control = scaled((unit**2*OFFSET/4-moment(control, h, 1))/pivot, even_mode)
    assert moment(a, h, 1) == moment(control, b, 1) == 0
    assert value(control, 1) == -2*unit and value(h, 1) == 8*unit
    assert 4*(moment(z_control, b, 1)+moment(control, h, 1)) == unit**2*OFFSET
    assert e[6] == Q(1, 6) and len(e) == 13 and e[-1] == -1/(12*d**2)
    assert all(len(poly) <= 11 for poly in (control, z_control))
    assert len(b) == 12 and len(h) <= 12
    assert value(e, d) > value(e, 1) > 0
    return a, e, b, h, control, z_control, unit, w


def falling_root(e, d, y):
    return root(lambda x: value(e, x)-y, d, 2*d, 140)


def increment(p, start, ends):
    q = primitive(p)
    return minus(poly_box(q, ends), (value(q, start),)*2)


def calibration_faces():
    sides = []
    for d in FACES:
        _, e, _, _, _, _, _, density = profile(d)
        raw = increment(density, d, falling_root(e, d, 0))
        sides.append(sign(raw))
        inside(raw, *({FACES[0]: ('-.000466', '-.000464'),
                       FACES[1]: ('.000777', '.000780')}[d]))
        for theta in (Q(-1, 10**5), Q(1, 10**5)):
            profile(d, theta)
    assert sides == [-1, 1]
    for d, expected in zip(TIGHT, (-1, 1)):
        a, e, b, _, _, _, _, raw = profile(d)
        assert sign(increment(raw, d, falling_root(e, d, 0))) == expected
        w = wronskian(a, b)
        assert sign(increment(w, d, falling_root(e, d, -Q(1, 2)))) < 0
    print("Opposite cap-action faces on two rational brackets; terminal orientation checked.")


def uniform_regular_turn():
    """Enclose the corrected turning comparison throughout the narrow d interval."""
    r = Q(1, 4)
    gamma = 1/TIGHT[1]**2, 1/TIGHT[0]**2
    weight = [0]*8 + multiply(SHOULDER, SHOULDER)
    numerator_mode = polynomial({0: 1, 2: 1})
    denominator_mode = polynomial({2: 3, 4: -1})
    numerators, denominators, roots = [], [], []
    for d in TIGHT:
        cap = polynomial({0: 1, 2: -1/d**2})
        omega = multiply(weight, multiply(cap, cap))
        numerators.append(integral(multiply(omega, numerator_mode), 1))
        denominators.append(integral(multiply(omega, denominator_mode), 1))
        _, e, *_ = profile(d)
        roots.append(falling_root(e, d, value(e, r)))
    # Both moments increase with d pointwise on (0,1). Also E(v)-E(r)
    # increases with d at fixed v>r, so its falling zero increases with d.
    k = -numerators[1]/denominators[0], -numerators[0]/denominators[1]
    v = roots[0][0], roots[1][1]

    def raw_primitive(z):
        values = []
        weights = [(ONE, ONE), times((-2, -2), gamma), times(gamma, gamma)]
        for mode in (numerator_mode, denominator_mode):
            base = multiply(weight, mode)
            result = (Q(0), Q(0))
            for i, coeff in enumerate(weights):
                h = primitive([0]*(2*i)+base)
                result = plus(result, times(coeff, poly_box(h, z)))
            values.append(result)
        return plus(values[0], times(k, values[1]))

    def slope(z):
        square = times(z, z)
        shoulder = minus((ONE, ONE), square)
        return times(times(times(square, square), times(shoulder, shoulder)),
                     minus((ONE, ONE), times(gamma, square)))

    falling = slope(v)
    assert falling[1] < 0
    rr = (r, r)
    # B^0 gives W=-raw_primitive; beta is positive.
    exterior = minus(raw_primitive(TIGHT), raw_primitive(v))
    quotient = times(exterior, (1/falling[1], 1/falling[0]))
    correction = minus(times((-1, -1), raw_primitive(rr)), times(slope(rr), quotient))
    inside(correction, '-3.1e-7', '-2.8e-7')
    print("At r=1/4 the corrected action divided by beta is in (-3.1e-7,-2.8e-7) throughout the narrow box.")


def finite_profile(d, theta, lam):
    a, _, b, h, control, z_control, unit, _ = profile(d, theta)
    eta = lam**2
    length = HIGHER[-2]/(6*HIGHER[-1]*b[5]*eta)
    height = 6*HIGHER[-1]*length**6
    even = primitive([0]+linear_combination([1, eta, lam**3], [a, control, z_control]))
    odd = linear_combination([eta, lam**3], [b, h])
    for i, coeff in enumerate(((LINEAR[0]+LINEAR[1])/2, *HIGHER[:3]), 1):
        (odd if i % 2 else even)[i] += coeff*length**i/height
    assert height*odd[5]/length**5 == HIGHER[-2]
    assert abs(height*even[6]/length**6-HIGHER[-1]) < Q(1, 10**16)
    assert len(even) == 13 and even[12] < 0 and len(odd) == 12 and odd[11]
    return even, odd, unit, 20/length


def finite_upper_actions(even, odd, s, d, lam, unit, cutoff):
    # Use the spatial variable v=(z-1)/lambda. Its limiting cubic has
    # the same unit as in (z^2/2-1/2)/lambda since dz/du=1 at z=1.
    f = linear_combination([1, s], [even, odd])
    dz = deriv(f)
    extent = Q(70)
    a = -2-s
    local = shift_at_one(f, lam)
    slope = scaled(1/(unit*lam**2), shift_at_one(dz, lam))
    remainder = difference(slope, [Q(a), Q(0), Q(3)])
    for p in (remainder, deriv(remainder)):
        assert sum(abs(c)*extent**i for i, c in enumerate(p)) < Q(1, 100)
    # The return reference uses only the rising part before the cap.
    for lo, hi in ((cutoff, Q(1, 4)), (Q(1, 4), Q(3, 4)),
                   (Q(3, 4), 1-extent*lam), (1+extent*lam, (1+d)/2)):
        assert bernstein_positive(dz, lo, hi)
    k = ONE if a == -3 else Q('.57735')
    def zero(p, target, left, right):
        return root(lambda x: value(p, x)-target, left, right)
    bottom_at = zero(slope, 0, k-Q('.02'), k+Q('.02'))
    peak_at = zero(slope, 0, -k-Q('.02'), -k+Q('.02'))
    bottom, peak = poly_box(local, bottom_at), poly_box(local, peak_at)
    assert 0 < bottom[0] <= bottom[1] < peak[0]
    action = primitive(multiply(dz, dz)[1:])
    action[0] -= value(action, cutoff)
    action = shift_at_one(action, lam)
    log_error = 3*dz[0]**2/cutoff
    landing_roots = [zero(local, t, -2*k-Q('.05'), -2*k+Q('.05')) for t in bottom]
    landing_action = poly_box(action, (landing_roots[0][0], landing_roots[1][1]))
    skipped = minus(poly_box(action, bottom_at), landing_action)
    assert skipped[0] > 0
    rows = []
    for sample in SAMPLES:
        y = value(even, 1)+unit*lam**3*sample
        assert all(y < lo or y > hi for lo, hi in (bottom, peak))
        if bottom[1] < y < peak[0]:
            result = landing_action
        else:
            ends = (-extent, peak_at[0]) if y < bottom[0] else (bottom_at[1], extent)
            entry = zero(local, y, *ends)
            result = poly_box(action, entry)
            if y > bottom[1]:
                result = minus(result, skipped)
        # The retained part of the discarded positive term dz[0]^2/z is bounded
        # by length times its maximum on [cutoff,3].
        rows.append((result[0], result[1]+log_error))
    return rows


def compact_returns(even, odd, cutoff, d):
    """Enclose true finite-field returns by sextic flows on the compact strip."""
    length = 20/cutoff
    height = 6*HIGHER[-1]*length**6
    mid = sum(LINEAR)/2
    seed = [0, mid, *HIGHER]
    seed_drift = Q(1, 10**60)
    full = linear_combination([1, 1], [even, odd])
    physical = [c*height/length**i for i, c in enumerate(full)]
    error = difference(physical, seed)
    assert error[0] == error[1] == 0
    assert sum(abs(c)*10**(i-1) for i, c in enumerate(error) if i) < seed_drift
    _, e, *_ = profile(d)
    reference = (value(e, 1)+value(e, d))/2
    wall = (1+9*d)/10
    assert reference*height > 9
    # The lower variable a also varies in the full shooting interval. Enclose
    # G_s(x)-8x^6 on the entire inward tail from the finite section L/2.
    for s in (1, -1):
        template = linear_combination([1, s], [even, odd])
        for a in LINEAR:
            excess = template[:]
            excess[1] += s*(a-mid)*length/height
            excess[6] -= 8*length**6/height
            for lo, hi in ((2/length, 10/length), (10/length, Q(1, 2))):
                assert bernstein_positive(excess, lo, hi)
            excess[6] += Q(79, 10)*length**6/height
            for lo, hi in ((Q(1, 2), ONE), (ONE, Q('1.1')), (Q('1.1'), wall)):
                assert bernstein_positive(excess, lo, hi)
            assert value(excess, wall)+wall**6*length**6/(10*height) > reference

    def endpoint_bound(a, orientation, xp, xm, yp, ym):
        # Levelwise ordering puts the two extreme fields at a +/- seed_drift.
        # orientation +1 encloses a lower bound on w_plus-w_minus.
        p = (Q(0), a+orientation*seed_drift, *HIGHER)
        reflected = tuple((-1)**i*c for i, c in enumerate(p))
        plus, _ = return_height(p, Q(xp), Q(yp))
        neg, _ = return_height(reflected, Q(xm), Q(ym))
        return plus[0]-neg[1] if orientation > 0 else neg[0]-plus[1]

    low, high = LINEAR
    for y, a, orient, xp, xm, margin in (
            ('.5', high, 1, 0, 0, '2e-7'), ('3.7', low, -1, 0, 0, '1e-24'),
            ('3.88', high, 1, 0, 0, '3e-27'), ('4.03', low, -1, 0, 0, '2e-28'),
            ('4.36', high, 1, 0, 0, '6e-32'), ('8', low, 1, 0, 2, '2e-32'),
            ('8', high, -1, 2, 0, '6e-32')):
        assert endpoint_bound(a, orient, xp, xm, y, y) > Q(margin)
    for y, orient, margin in ((Q(4, 5), -1, '5e-8'), (Q(0), 1, '.003'),
                              (Q(-20), -1, '.001')):
        # At x=2 each ordinate lies between its section value minus TAIL
        # and that value. Lower magnitude decreases in this inward ordinate.
        yp = y if orient > 0 else y-TAIL
        ym = y-SHIFT-TAIL if orient > 0 else y-SHIFT
        assert endpoint_bound(high if orient > 0 else low, orient, 2, 2, yp, ym) > Q(margin)
    print("Positive-tilt comparison profile: five compact signs, two shooting faces, three inward section comparisons.")


def shoulder_and_ledger():
    for i, (height, band) in enumerate(zip(SAMPLES, WIDTHS)):
        action, _ = outward_signs(height, ONE)
        inside(action, *band)
        assert sign(action) == (-1)**i
    d = Q(6, 5)
    lam = Q(1, 10**40)
    even, odd, unit, cutoff = finite_profile(d, Q(1, 10**20), lam)
    compact_returns(even, odd, cutoff, d)
    plus, neg = [finite_upper_actions(even, odd, s, d, lam, unit, cutoff) for s in (1, -1)]
    for i, (p, n, band) in enumerate(zip(plus, neg, WIDTHS)):
        diagnostic = times((-1/(unit**2*lam**5),)*2, minus(p, n))
        inside(diagnostic, *band)
        assert sign(minus(p, n)) == -(-1)**i
    compact = [1, -1, 1, -1, 1]
    shoulder = [sign(minus(p, n)) for p, n in zip(plus, neg)]
    signs = compact+shoulder+[-1, 1]+[-1, 1, -1]+[1, -1]
    assert len(signs) == 20 and all(p == -n for p, n in zip(signs, signs[1:]))
    assert len(signs)-1 == 19 > 18
    print("Full rational degree-12 profile passes both folds, landings and eight retained upper actions.")
    print("Conditional count using the analytic matches and detuning: 20 alternating signs, 19 brackets.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run Python with assertions enabled.")
    calibration_faces()
    uniform_regular_turn()
    shoulder_and_ledger()
