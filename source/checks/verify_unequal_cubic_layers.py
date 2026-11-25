"""Rational tests for the sextic graft, unequal cubic weights and degree cost."""

from fractions import Fraction as Q
import sys

from verify_outer import value
from verify_remote_odd_powers import bernstein_positive, polynomial, shifted
from verify_seed import HIGHER, LINEAR, multiply
from verify_square_shoulders import (
    SAMPLES, WIDTHS, cubic_rows, deriv, inside, linear_combination, minus,
    moment, monomial, poly_box, primitive, root, shift_at_one, solve, times,
)


OFFSET = Q(493, 2)
RHO = 1 + Q(1, 10**9)
FAR = Q(10**21)
INCOMING = (0, 3, 6, 7)
FACTOR = polynomial({0: 1, 2: -1})


def sign(box):
    assert box[0] * box[1] > 0, box
    return 1 if box[0] > 0 else -1


def scaled(c, p):
    return [c * x for x in p]


def outward_signs(height, rho):
    pos = cubic_rows(-1, Q(-8), height)
    neg = cubic_rows(-3, Q(8), height)
    return [minus(minus(times((rho, rho), p), n), (OFFSET, OFFSET))
            for p, n in zip(pos, neg)]


def cubic_limits():
    for weight in (Q(1), RHO):
        for i, (height, band) in enumerate(zip(SAMPLES, WIDTHS)):
            upper, inc = outward_signs(height, weight)
            if weight == 1:
                inside(upper, *band)
            assert sign(upper) == (-1)**i
            if i in INCOMING:
                assert sign(inc) == (-1)**INCOMING.index(i)
    for height, expected in ((-FAR, -1), (FAR, 1)):
        assert [sign(box) for box in outward_signs(height, RHO)] == [expected]*2


def seeded_profile(lam):
    eta = lam**2
    k = Q(-45, 4)
    r = scaled(Q(3, 2), multiply(monomial(2), FACTOR))
    r2 = multiply(r, r)
    b = scaled(k, multiply(multiply(monomial(5), FACTOR),
                           polynomial({0: 1, 2: Q(-17, 15)})))
    h = scaled(Q(1, 5), polynomial({7: -594, 9: 714}))
    d = scaled(Q(1, 13), polynomial({4: 342, 6: -420}))
    control = scaled(Q(9578943, 2912), monomial(4))
    assert all(moment(f, g, 1) == 0 for f, g in ((r2, b), (r2, h), (d, b)))
    assert 4*(moment(control, b, 1) + moment(d, h, 1)) == 9*OFFSET
    assert [value(deriv(g), 1) for g in (r, b)] == [-3, -3]
    assert value(d, 1) == -6 and value(h, 1) == 24
    raw = linear_combination([1, -1], [multiply(r2, deriv(b)), multiply(deriv(r2), b)])
    expected = scaled(Q(9, 16), multiply(monomial(9), multiply(
        multiply(FACTOR, multiply(FACTOR, FACTOR)), polynomial({0: -5, 2: 3}))))
    assert primitive(raw) == expected and value(expected, 2) == -54432
    cumulative = primitive([0]+b)
    assert cumulative == scaled(k, polynomial({7: Q(1,7), 9: Q(-32,135), 11: Q(17,165)}))
    assert Q(32,135)**2 < 4*Q(1,7)*Q(17,165)
    e = primitive([0] + r2)
    assert e[6] == Q(3, 8)
    assert [e[i]/e[6] for i in (6, 8, 10)] == [1, Q(-3, 2), Q(3, 5)]
    assert 1-Q(3, 2)**2/(4*Q(3, 5)) == Q(1, 16) and HIGHER[-1]/16 > 1
    assert b[5] == k and b[-1] < 0 and HIGHER[-2] < 0
    slope = linear_combination([1, eta, lam**3], [r2, d, control])
    even = primitive([0] + slope)
    odd = linear_combination([eta, lam**3], [b, h])
    assert len(even) == 11 and len(odd) == 10 and even[-1] > 0 > odd[-1]
    return even, odd


def weighted_profile(lam):
    eta = lam**2
    delta = 3*(1-RHO**3)/(1+RHO**3)
    assert -Q(1, 100) < delta < 0
    unit = 24+8*delta
    assert 24-8*delta == unit*RHO**3
    square = multiply(FACTOR, FACTOR)
    base = multiply(square, polynomial({0: 1, 4: 17}))
    odd0 = scaled(delta, multiply(monomial(3), multiply(square, FACTOR)))
    assert moment(base, odd0, 1) == 0
    assert moment(base, multiply(monomial(3), square), 1) == Q(128, 4095)
    slope = list(base)
    slope[0] += eta*unit*(-RHO-3)/2
    modes = [monomial(i) for i in (3, 5, 7)]
    rows = [[value(p, 1) for p in modes], [value(deriv(p), 1) for p in modes],
            [4*moment(slope, p, 1) for p in modes]]
    rhs = [-8*unit*lam**3, eta*unit*(3-RHO)/2,
           -unit**2*lam**5*OFFSET-4*moment(slope, odd0, 1)]
    weights = solve(rows, rhs)
    odd = linear_combination([1, *weights], [odd0, *modes])
    assert value(odd, 1) == rhs[0] and value(deriv(odd), 1) == rhs[1]
    assert 4*moment(slope, odd, 1) == -unit**2*lam**5*OFFSET
    even = primitive([0] + slope)
    assert even[2] > 0 and len(even) == 11 and even[-1] > 0
    assert odd[3] < 0 and len(odd) == 10 and odd[-1] == -delta
    return even, odd, unit


class FiniteLayer:
    """Two folds, one landing, and retained actions in v=(z-1)/lambda."""

    def __init__(self, even, odd, side, lam, unit, alpha, cubic_a, extent, cutoff):
        self.extent = extent
        self.height_unit = unit*lam**3
        self.center = value(even, 1)
        profile = linear_combination([1, side], [even, odd])
        self.global_profile = profile
        dz = deriv(profile)
        self.profile = shift_at_one(profile, lam)
        primitive_action = primitive(multiply(dz, dz)[1:])
        primitive_action[0] -= value(primitive_action, cutoff)
        self.global_action = primitive_action
        self.action = shift_at_one(primitive_action, lam)
        # The discarded integrand is dz[0]^2/z. On [cutoff,3] its
        # contribution to any retained action lies in [0,3*dz[0]^2/cutoff].
        self.log_error = 3*dz[0]**2/cutoff if cutoff else Q(0)
        assert cutoff or dz[0] == 0
        slope = scaled(1/(unit*lam**2), shift_at_one(dz, lam))
        err = linear_combination([1, -1], [slope, [alpha*cubic_a, 0, 3*alpha**3]])
        for p in (err, deriv(err)):
            assert sum(abs(c)*extent**i for i, c in enumerate(p)) < Q(1, 100)
        # Certify the positive radial derivative on the full retained ray.
        radial = dz[:]
        while radial[0] == 0:
            radial = radial[1:]
        for low, high in ((cutoff, Q(1, 4)), (Q(1, 4), Q(3, 4)),
                          (Q(3, 4), 1-extent*lam)):
            assert bernstein_positive(radial, low, high)
        assert all(c > 0 for c in shifted(radial, 1+extent*lam))
        k = (Q(1) if cubic_a == -3 else Q(".57735"))/alpha

        def zero(poly, target, lo, hi):
            return root(lambda x: value(poly, x)-target, lo, hi)

        self.minimum = zero(slope, 0, k-Q(".02"), k+Q(".02"))
        self.maximum = zero(slope, 0, -k-Q(".02"), -k+Q(".02"))
        self.bottom = poly_box(self.profile, self.minimum)
        self.peak = poly_box(self.profile, self.maximum)
        assert 0 < self.bottom[0] and self.bottom[1] < self.peak[0]
        lo, hi = -2*k-Q(".05"), -2*k+Q(".05")
        left, right = (zero(self.profile, y, lo, hi) for y in self.bottom)
        self.landing = poly_box(self.action, (left[0], right[1]))
        self.skip = minus(poly_box(self.action, self.minimum), self.landing)
        assert self.skip[0] > 0

    def at(self, sample):
        height = self.center+self.height_unit*sample
        assert all(height < box[0] or height > box[1] for box in (self.bottom, self.peak))
        lo, hi = ((-self.extent, self.maximum[0]) if height < self.bottom[0]
                  else (self.minimum[1], self.extent))
        entry = root(lambda x: value(self.profile, x)-height, lo, hi)
        incoming = poly_box(self.action, entry)
        if height > self.bottom[1]:
            incoming = minus(incoming, self.skip)
        upper = self.landing if self.bottom[1] < height < self.peak[0] else incoming
        return [(lo, hi+self.log_error) for lo, hi in (upper, incoming)]

    def beyond_folds(self, height):
        assert height > self.peak[1]
        entry = root(lambda z: value(self.global_profile, z)-height, Q("1.9"), Q("2.1"), 170)
        lo, hi = minus(poly_box(self.global_action, entry), self.skip)
        return lo, hi+self.log_error


def finite_samples():
    for weighted in (False, True):
        lam = Q(1, 10**45) if weighted else Q(1, 10**12)
        if weighted:
            even, odd, unit = weighted_profile(lam)
            settings = ((1, RHO, -1), (-1, Q(1), -3))
            extent = Q(2*10**7)
            cutoff = Q(0)
        else:
            even, odd = seeded_profile(lam)
            scale_ratio = Q(8,3)*HIGHER[-1]
            length = HIGHER[-2]/(Q(-45,4)*scale_ratio*lam**2)
            height = scale_ratio*length**6
            for power, coeff in enumerate(((LINEAR[0]+LINEAR[1])/2, *HIGHER[:3]), 1):
                (odd if power % 2 else even)[power] += coeff*length**power/height
            cutoff = 20/length
            unit = Q(3)
            settings = ((1, Q(1), -3), (-1, Q(1), -1))
            extent = Q(70)
        sides = [FiniteLayer(even, odd, s, lam, unit, alpha, a, extent, cutoff)
                 for s, alpha, a in settings]
        orientation = 1 if weighted else -1
        for i, (t, band) in enumerate(zip(SAMPLES, WIDTHS)):
            left, right = (s.at(t) for s in sides)
            scale = orientation/(unit**2*lam**5)
            u, inc = [times((scale, scale), minus(p, n)) for p, n in zip(left, right)]
            if not weighted:
                inside(u, *band)
            assert sign(u) == (-1)**i
            if i in INCOMING:
                assert sign(inc) == (-1)**INCOMING.index(i)
        if weighted:
            for t, expected in ((-FAR, -1), (FAR, 1)):
                left, right = (s.at(t) for s in sides)
                assert all(sign(minus(p, n)) == expected for p, n in zip(left, right))
        else:
            r = scaled(Q(3, 2), multiply(monomial(2), FACTOR))
            height = value(primitive([0]+multiply(r, r)), 2)
            left, right = (s.beyond_folds(height) for s in sides)
            inside(times((1/(2*lam**2),)*2, minus(left, right)), -54433, -54431)
    print("Both finite degree-ten profiles pass two-fold, landing and retained-action tests.")


def weighted_inward_flight():
    # z=R*v, entrance ordinate -S*R^10; allow every 0<=a<=1 and L>=L0.
    even, odd, _ = weighted_profile(Q(1, 10**45))
    radius, length = Q(10**8), Q(10**200)
    top, e2 = odd[9], even[2]
    assert 1 < even[10] < 2
    # The normalized denominators h_s=(R^10+G_s(R*v)+s*a*e2*R*v/L)/R^10
    # lie between (1+v^10)/2 and 3*(1+v^10) for the entire ray.
    lower_terms = sum((abs(x)+abs(y))*radius**(i-10)
                      for i, (x, y) in enumerate(zip(even[:10], odd)))
    linear = e2*radius**-9/length
    assert lower_terms+linear < Q(1, 2)
    # The scaled odd numerator divided by top*R^9 is v^9 plus an error
    # bounded by odd_error*(1+v^9), including the shooting parameter.
    odd_error = sum(abs(c)*radius**(i-9)/top for i, c in enumerate(odd[:9]))
    odd_error += e2/(length*top*radius**8)
    margin = 1/Q(9*(1+2**10)**2)-8*odd_error
    # Restrict the positive leading integral to [1,2]. Over the full ray,
    # integral v*(1+v^9)/(h_+*h_-) is <8 by splitting at 1.
    assert margin > 0
    slow = e2**2/length**2*radius**-18
    # Each inward increment in rescaled ordinate is at most 2*slow.
    # Replacing h_s+increment by h_s in its integral costs <5*slow per side.
    assert 10*slow < 2*top*margin/radius
    print("The full negative incoming difference is strictly negative for L>=10^200, uniformly in shooting.")


def support_obstruction():
    # Verify all coefficients of the height and prefix constraints by polarization.
    def integral(p):
        return sum(2*c/Q(i+1) for i, c in enumerate(p) if i % 2 == 0)

    square = multiply(FACTOR, FACTOR)
    fourth = multiply(square, square)
    for e in (0, 1):  # affine in e
        tpe = [e, Q(1)]
        for a, b, c in ((1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1)):
            w = [a, b, c]
            assert integral(multiply(tpe, multiply(square, w))) == Q(16,105)*(7*a*e+b+c*e)
            target = 143*a*a*e+26*a*b+26*a*c*e+13*b*b*e+6*b*c+3*c*c*e
            assert integral(multiply(tpe, multiply(fourth, multiply(w,w)))) == Q(256,45045)*target
    assert 10**2-4*3*39 < 0


def ledger():
    def changes(xs):
        return sum(a != b for a, b in zip(xs, xs[1:]))

    finite = [1,-1,1,-1,1]
    layer = [-1,1]*4
    assert changes([*finite, *layer, -1, 1]) == 14
    assert changes([1, -1, 1, -1, 1, -1, 1]) == 6  # remote to physical 1, 0
    assert max(10+2, 9+4) == 13 and 14+6+2 == 22
    weighted = [-1, 1]*5
    assert changes([*weighted, -1]) == 10
    assert changes([*([-1, 1]*3), -1]) == 6
    for t in (0, 1, 2, 120):
        m, j = 10+2*t, 9+2*t
        assert max(m+2, j+4) == 13+2*t and max(m+4, j+4) == 14+2*t
        q, r = 10+3*t, 6+4*t
        assert q+r+2 == 18+7*t == 15+7*(1+t)-4
        assert (q+1)+(r+2)+2 == 21+7*t == 18+7*(1+t)-4
    print("Seeded core: 14 at degree 10 and 22 at 13. Weighted core: 10 at 10 and 18 at 13.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run with assertions enabled.")
    cubic_limits()
    finite_samples()
    weighted_inward_flight()
    support_obstruction()
    ledger()
