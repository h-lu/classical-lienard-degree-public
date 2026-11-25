"""Rational checks for the quadratic core, far-return signs and exact degree ledger."""

from fractions import Fraction as Q
from math import comb
import sys

from verify_seed import multiply
from verify_outer import value
from verify_square_shoulders import (
    SAMPLES, WIDTHS, deriv, inside, linear_combination, minus, moment, monomial,
    plus, poly_box, primitive, root, shift_at_one, solve, times,
)


def polynomial(terms):
    coefficients = [Q(0)]*(max(terms)+1)
    for power, coefficient in terms.items():
        coefficients[power] = Q(coefficient)
    return coefficients


def shifted(poly, center):
    return [sum(coefficient*comb(i,k)*center**(i-k)
                for i,coefficient in enumerate(poly) if i >= k)
            for k in range(len(poly))]


def bernstein_positive(poly, a, b):
    """Prove positivity on [a,b] using shifted Bernstein coefficients."""
    coeff = [c*(b-a)**i for i,c in enumerate(shifted(poly,a))]
    n = len(coeff)-1
    return all(sum(coeff[i]*Q(comb(k,i),comb(n,i)) for i in range(k+1)) > 0
               for k in range(n+1))


def scale(poly, radius, divisor):
    return [c*radius**i/divisor for i,c in enumerate(poly)]


def example():
    lam = Q(1,10**12)
    eta = lam**2
    b = polynomial({3:"1.5", 5:"-1.5"})
    h = polynomial({5:-1, 7:Q(85,77)})
    d = polynomial({0:2, 8:88, 10:-96})
    cap = Q(619003,16)
    zcontrol = polynomial({8:-cap, 10:cap})
    g = 231
    assert moment(multiply(b,b),h,1) == 0
    assert g*value(h,1) == 24 and value(d,1) == -6
    assert moment(d,b,1) == 0
    assert 4*(moment(zcontrol,b,1)+g*moment(d,h,1)) == 9*Q(493,2)
    slope = linear_combination([1,eta,eta*lam], [multiply(b,b),d,zcontrol])
    even = primitive([Q(0)]+slope)
    odd = linear_combination([eta,g*eta*lam], [b,h])
    b7, a12 = odd[7], even[12]
    assert even[2] == eta and a12 > 0 and b7 > 0
    # The four added terms occur in the order odd, odd, even, odd.
    b9 = -b7/Q(10**40)
    b11 = b9/Q(10**80)
    a14 = a12/Q(10**120)
    b13 = b11/Q(10**160)
    even = linear_combination([1,1], [even,polynomial({14:a14})])
    odd = linear_combination([1,1], [odd,polynomial({9:b9, 11:b11, 13:b13})])
    return lam, even, odd


def two_inflections():
    roots = (Q(1), Q(16))
    b = polynomial({3:Q(3,2)})
    for q in roots:
        b = multiply(b, [Q(1), Q(0), -1/q**2])
    square = multiply(b,b)
    modes = [monomial(i) for i in (5,7,9)]
    rows = [[moment(square,p,q) for p in modes] for q in roots]
    odd = linear_combination([-1,*solve([r[1:] for r in rows],[r[0] for r in rows])],modes)
    assert all(moment(square,odd,q) == 0 for q in roots)
    assert all(value(odd,q)*value(deriv(b),q) < 0 for q in roots)
    assert odd[-1]*b[-1] < 0
    modes = [monomial(i) for i in (8,10,12,14)]
    rows = [[value(p,q) for p in modes] for q in roots]
    rows += [[moment(p,b,q) for p in modes] for q in roots]
    coeffs = solve(rows,[Q(-8),Q(0),Q(0),Q(0)])
    d = linear_combination([1,*coeffs],[[Q(2)],*modes])
    assert [value(d,q) for q in roots] == [-6,2]
    assert all(moment(d,b,q) == 0 for q in roots)
    # Two independent offset targets including the full D H' correction.
    targets = [Q(1),Q(-3)]
    rhs = [target-moment(d,odd,q) for target,q in zip(targets,roots)]
    zcontrol = linear_combination(solve(rows,[Q(0),Q(0),*rhs]),modes)
    assert all(value(zcontrol, q) == 0 for q in roots)
    assert [moment(zcontrol,b,q)+moment(d,odd,q) for q in roots] == targets
    print("Two separated cubic-germ inflections satisfy the exact odd moments and affine even constraints.")


class Side:
    """Exact critical, entry and landing enclosures for the one-shoulder profile."""
    def __init__(self, sign, lam, even, odd):
        self.lam = lam
        self.height_unit = 3*lam**3
        self.center = value(even,1)
        self.profile = linear_combination([1,sign],[even,odd])
        derivative = deriv(self.profile)
        self.action = primitive(multiply(derivative,derivative)[1:])
        self.local_profile = shift_at_one(self.profile,lam)
        self.local_action = shift_at_one(self.action,lam)
        local_slope = shift_at_one(derivative,lam)
        a = -2-sign
        residual = linear_combination([1/(3*lam**2),-1],[local_slope,[a,0,3]])
        for p in (residual,deriv(residual)):
            assert sum(abs(c)*70**i for i,c in enumerate(p)) < Q(1,100)
        # Outside the two normalized folds the physical slope is positive.
        slope_z = derivative[1:]  # G'(z)/z
        for left,right in ((Q(0),Q(1,4)),(Q(1,4),1-70*lam)):
            assert bernstein_positive(slope_z,left,right)
        assert all(c > 0 for c in shifted(slope_z,1+70*lam))
        k = Q(1) if a == -3 else Q(".57735")
        self.min_root = root(lambda x:value(local_slope,x),k-Q(".02"),k+Q(".02"))
        self.max_root = root(lambda x:value(local_slope,x),-k-Q(".02"),-k+Q(".02"))
        self.bottom = poly_box(self.local_profile,self.min_root)
        self.peak = poly_box(self.local_profile,self.max_root)
        assert self.bottom[0] > 0 and self.bottom[1] < self.peak[0]
        left, right = -2*k-Q(".03"), -2*k+Q(".03")
        land_lo = root(lambda x:value(self.local_profile,x)-self.bottom[0],left,right)
        land_hi = root(lambda x:value(self.local_profile,x)-self.bottom[1],left,right)
        self.landing = poly_box(self.local_action,(land_lo[0],land_hi[1]))
        self.skip = minus(poly_box(self.local_action,self.min_root),self.landing)
        assert self.skip[0] > 0

    def local_actions(self, height):
        target = self.center+self.height_unit*height
        for window_edge in (self.bottom,self.peak):
            assert target < window_edge[0] or target > window_edge[1]
        if target < self.bottom[0]:
            low, high = Q(-70), self.max_root[0]
        else:
            low, high = self.min_root[1], Q(70)
        entry = root(lambda x:value(self.local_profile,x)-target,low,high)
        incoming = poly_box(self.local_action,entry)
        if target > self.bottom[1]:
            incoming = minus(incoming,self.skip)
        upper = self.landing if self.bottom[1] < target < self.peak[0] else incoming
        return upper, incoming

    def remote_action(self, radius, height, relative):
        low, high = radius*(1-2*relative), radius*(1+2*relative)
        assert low > 2 and height > self.peak[1]
        bracket = root(lambda x:value(self.profile,x)-height,low,high,32)
        return minus(tuple(value(self.action,z) for z in bracket),self.skip)


def action_signs(lam, even, odd, sides, radii):
    for index,(sample,band) in enumerate(zip(SAMPLES,WIDTHS)):
        pos, neg = [side.local_actions(sample) for side in sides]
        upper = times((1/(9*lam**5),)*2,minus(neg[0],pos[0]))
        inside(upper,*band)
        if index in (0,3,6,7):
            desired = 1 if index in (0,6) else -1
            assert min(desired*x for x in minus(neg[1],pos[1])) > 0
    # At regular levels beyond the folds both actions retain the same branches.
    signs = (-1,1,-1,1,-1)
    leading_ratios = (Q(-72,17), Q(-24,19), Q(8,7), Q(-28,23), Q(28,25))
    for radius,desired,ratio in zip(radii,signs,leading_ratios):
        height = value(even,radius)
        odd_height = value(odd,radius)
        relative = abs(odd_height)/height
        assert 0 < relative < Q(1,100)
        p, n = [side.remote_action(radius,height,relative) for side in sides]
        difference = minus(p,n)
        assert min(desired*x for x in difference) > 0
        normalized = times(difference,(radius**2/(height*odd_height),)*2)
        inside(normalized,ratio-Q(1,10**6),ratio+Q(1,10**6))
    print("One exact profile has all moving-height signs and five distant action signs -,+,-,+,-.")


def integral_margin(even, odd, sides, radius, expected):
    """Bound integral of z*C / ((w+G_+)(w+G_-)) on the entire positive ray."""
    w = value(even,radius)
    magnitude = sum(abs(c)*radius**i for i,c in enumerate(odd))
    numerator = [Q(0)]+scale(odd,radius,magnitude)
    side_polys = [scale(side.profile,radius,w) for side in sides]
    end = 16
    assert radius*end >= 2
    tail_bounds = []
    for degree, leading in ((12,Q(1,20)),(14,even[14]/2)):
        inverse = w/(leading*radius**degree)
        tail_bounds.append(inverse**2*sum(
            abs(c)*Q(end)**(k+2-2*degree)/(2*degree-k-2)
            for k,c in enumerate(numerator[1:]) if c))
    tail = min(tail_bounds)
    # Split [0,16] on geometric panels. The fields are positive, so the
    # normalized denominators have the additional lower bound 1 on every panel.
    total = (Q(0),Q(0))
    grid = 2**100
    edges = [Q(0),Q(1,64)]
    while edges[-1] < end:
        start = edges[-1]
        edges.extend(start*(1+Q(i,48)) for i in range(1,49))
    assert edges[-1] == end
    for left,right in zip(edges,edges[1:]):
        iv = (left,right)
        low_denom = high_denom = Q(1)
        for poly in side_polys:
            low, high = poly_box(poly,iv)
            low_denom *= max(Q(1),1+low)
            high_denom *= 1+high
        interval = times(poly_box(numerator,iv),
                         ((right-left)/high_denom,(right-left)/low_denom))
        lower = (interval[0]*grid).__floor__()
        upper = (interval[1]*grid).__ceil__()
        total = plus(total,(Q(lower,grid),Q(upper,grid)))
    margin = min(expected*t for t in total)-tail
    assert margin > Q(1,20)
    return w, margin*radius**2*magnitude/w**2


def physical_negative(lam, even, odd, sides, far):
    # The polynomial residual, expanded about 2, is strictly positive.
    lower_tail = polynomial({12:Q(1,20),14:even[14]/2})
    for side in sides:
        residual = linear_combination([1,-1],[side.profile,lower_tail])
        assert all(v > 0 for v in shifted(residual,Q(2)))
    radii = (Q(1,2),Q(3),*far)
    kernel_signs = (1,-1,1,-1,-1,-1,-1)
    eta, space = lam**2, 10**3000
    eps = (eta/space)**2
    for radius, expected in zip(radii,kernel_signs):
        w, margin = integral_margin(even,odd,sides,radius,expected)
        lower = min(w/(1+2**14), even[14]/4)
        relative = (eta/space)/lower + 2*eps/lower**2
        assert relative < Q(1,2) and 4*eps/lower**2 < 1
        # Absolute error in (I_+-I_-)/(S*eps), compared with minus twice
        # the signed integral. The perturbation includes every a in [-1,0].
        assert 4*relative/lower < 2*margin
    print("Seven true negative entrance comparisons at L=10^3000 have signs -,+,-,+,+,+,+.")


def shooting(lam, even, odd, sides):
    eta, space = lam**2, Q(10**3000)
    for side in sides:
        residual = list(side.profile[2:])  # G_s(z)/z^2 - eta/2
        residual[0] -= eta/2
        for low,high in ((Q(0),Q(1,4)),(Q(1,4),Q(1)),(Q(1),Q(2))):
            assert bernstein_positive(residual,low,high)
        assert all(c > 0 for c in shifted(residual,Q(2)))
    # F(sx) >= x^2/2-x >= v/2-1, so a reference arriving at positive
    # y satisfies v=x^2/2 <= 2y+6. Its linear memory on [1,Y] is <=exp(1-Y).
    height, xmax = 16000, Q(255)
    assert 2*(2*height+6) < xmax**2
    even_error = sum(i*abs(c)*(xmax/space)**(i-2)/eta
                     for i,c in enumerate(even) if i > 2)
    odd_error = sum(abs(c)*(xmax/space)**(i-3)/eta
                    for i,c in enumerate(odd) if i > 3)
    assert even[2] == eta and even_error < 1
    assert odd[3]/eta == Q(3,2) and odd_error < Q(1,2)
    # Thus the physical common v-slope is in [1,3], and the odd term
    # excluding ax lies between x^3/L and 2x^3/L up to xmax.
    small_x = Q(1,4)
    assert Q(3,2)*small_x**2 + small_x + 2*small_x**3/space < Q(1,2)
    assert 2*xmax**2/space < Q(1,2)
    # On y in [1,2] the lower x-barrier is 1/4. Even using e<3, e>2,
    # the integrated odd forcing at a=0 dominates the incoming memory.
    memory = Q(2*height+6, 2**(height-1))
    assert memory < 2*small_x**3/(27*space) and memory < small_x/27
    print("At the same L the physical reference has opposite shooting signs at a=-1 and a=0.")


def changes(signs):
    return sum(x != y for x,y in zip(signs,signs[1:]))


def ledger(h, e):
    m = 8+4*h
    j = 5+2*h
    upper, pos = [1], []
    for i in range(1,h+1):
        s = (-1)**i
        upper.extend(s*(-1)**t for t in range(8))
        pos.extend(-s*(-1)**t for t in range(4))
    neg = [-(-1)**i for i in range(h+1)]
    b = -(-1)**h
    for step in range(h+1):
        if step:
            j += 2
            b = -b
        assert j <= m-3
        upper.append(-b)
        pos.append(b)
        neg.append(-b)
    assert j == m-3 and b == -1
    for step in range(e):
        if step == 1:
            m += 2
        else:
            j += 2
        far_upper = b if j == m-1 else -b
        upper.append(far_upper)
        pos.append(-far_upper)
    upper.append(-upper[-1])
    entrance = pos[::-1]+neg
    if m-j == 3:
        entrance.append(-1)
    q, r = changes(upper),changes(entrance)
    assert q == 8*h+3+e and r == 6*h+(4,4,6,6)[e]
    direct = b if m-j == 1 else -b
    count = q+r+1+int(entrance[-1] != direct)
    n = max(m+2,j+4)
    assert n == 10+4*h+e and count == 14*h+(9,10,13,14)[e]
    return n,count


def degree_comparison():
    for h in range(1, 150):
        for e in range(4):
            n, count = ledger(h, e)
            assert 2*count >= 7*n - 63
    print("The oriented lists and tail give the four exact-degree counts of Proposition 5.1 for 1 <= h <= 149.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run with assertions enabled.")
    two_inflections()
    lam, even, odd = example()
    sides = [Side(s,lam,even,odd) for s in (1,-1)]
    far = tuple(Q(10**k) for k in (10,30,50,70,90))
    action_signs(lam,even,odd,sides,far)
    physical_negative(lam,even,odd,sides,far)
    shooting(lam,even,odd,sides)
    degree_comparison()
