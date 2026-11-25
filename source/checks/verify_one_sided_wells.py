"""Exact rational certificates for the octic core and two successive one-sided wells."""

from fractions import Fraction as Q
from math import comb
import sys

from verify_outer import value
from verify_seed import multiply
from verify_square_shoulders import (
    SAMPLES, WIDTHS, deriv, inside, linear_combination, minus, moment, plus,
    poly_box, primitive, root, shift_at_one, times,
)
from verify_remote_odd_powers import changes, polynomial, scale, shifted


def sign_interval(interval, sign):
    assert min(sign*x for x in interval) > 0


def balanced_core(lam):
    eta = lam**2
    r = polynomial({1:Q(3,2), 3:Q(-3,2)})
    b = polynomial({3:5, 5:Q(-23,2), 7:Q(13,2)})
    h = polynomial({5:1, 7:Q(-65,49)})
    d = polynomial({0:2, 2:Q(304,17), 4:Q(-440,17)})
    cap, g = Q(60166953,4624), Q(147,2)
    z = polynomial({2:cap, 4:-cap})
    r2 = multiply(r,r)
    for f, odd in ((r2,b),(r2,h),(d,b)):
        assert moment(f,odd,1) == 0
    assert 4*(moment(z,b,1)+g*moment(d,h,1)) == -9*Q(493,2)
    assert [value(deriv(r),1),value(deriv(b),1),g*value(h,1),value(d,1)] == [-3,3,-24,-6]
    wronskian = primitive(linear_combination(
        [1,-1], [multiply(r2,deriv(b)),multiply(deriv(r2),b)]))
    product = polynomial({5:Q(9,8)})
    for constant in (-1,-1,-1,-2):
        product = multiply(product,polynomial({0:constant,2:1}))
    assert product == wronskian
    assert value(wronskian, 2) == 1944
    even = primitive([Q(0)]+linear_combination([1,eta,eta*lam],[r2,d,z]))
    odd = linear_combination([eta,g*eta*lam],[b,h])
    assert even[2] == eta and odd[3] == 5*eta and odd[7] > 0
    assert (len(even)-1,len(odd)-1) == (8,7)
    return even,odd


def action_primitive(poly):
    slope = deriv(poly)
    return primitive(multiply(slope,slope)[1:])


def endpoint_table():
    bands = [(22,23),(-9,-8),(-36,-35),(-61,-60),(-82,-81),(-102,-101),(-119,-118)]
    for d,band in enumerate(bands,3):
        m = 2*d
        t = Q(d,d+1)**d/(d+1)
        inner = root(lambda x:x**d*(1+x)-t,Q(0),Q(1),64)
        outer = root(lambda x:x**d*(x-1)-t,Q(1),Q(2),64)
        primitives = [action_primitive(polynomial({m:1,m+1:2*s,m+2:1})) for s in (1,-1)]
        p = tuple(value(primitives[0],r) for r in inner)
        n = tuple(value(primitives[1],r)-value(primitives[1],1) for r in outer)
        inside((100*(1-n[1]/p[0]),100*(1-n[0]/p[1])),*band)
        print(f"well endpoint m={m}: 100 B_{d} in ({band[0]}, {band[1]})")
    assert 5*11*Q(2,3)**10 < 1 and Q(2,3)*12/11 < 1
    assert 2+Q(9,20) > Q(9,4)
    print("The outer-action enclosures are negative for m=8 through 18 and positive for m=6.")


def octic_well_difference(a, level):
    """Enclose the retained action difference at a level above the peak of the minus profile."""
    up, down = [polynomial({8: 1, 9: side*a, 10: 1}) for side in (1, -1)]

    def critical(x):
        return 8 - 9*a*x + 10*x*x

    peak = root(critical, Q('0.7'), Q('0.9'))
    bottom = root(critical, Q('0.9'), Q(1))
    assert poly_box(down, peak)[1] < level
    min_level = poly_box(down, bottom)
    assert min_level[0] > 0
    landing_lo, _ = root(lambda x: value(down, x)-min_level[0], Q(0), peak[0])
    _, landing_hi = root(lambda x: value(down, x)-min_level[1], Q(0), peak[0])
    inner = root(lambda x: value(up, x)-level, Q(0), Q(1))
    outer = root(lambda x: value(down, x)-level, Q(1), Q(2))
    up_action, down_action = map(action_primitive, (up, down))
    retained = plus(poly_box(down_action, (landing_lo, landing_hi)),
                    minus(poly_box(down_action, outer), poly_box(down_action, bottom)))
    return minus(poly_box(up_action, inner), retained)


def figure_action():
    level = Q('0.0076')
    inside(octic_well_difference(Q('1.998'), level), '0.00038127', '0.00038129')
    inside(octic_well_difference(Q('1.9998'), level), '-0.000017', '-0.000016')
    print("Nearby well at a=1.998, Y=0.0076: positive retained action difference in (0.00038127,0.00038129).")
    print("Figure 2 at a=1.9998, Y=0.0076: the retained action difference is in (-0.000017,-0.000016).")


def variations(coeffs):
    nonzero = [x > 0 for x in coeffs if x]
    return changes(nonzero)


def shape_certificate(poly, cuts, root_count):
    """Descartes on each bounded interval by Bernstein coefficients, then on its tail."""
    count = 0
    n = len(poly)-1
    endpoint_values = [value(poly,x) for x in cuts]
    assert all(endpoint_values)
    if root_count == 0:
        assert min(endpoint_values) > 0
    for a,b in zip(cuts,cuts[1:]):
        weights = [c*(b-a)**i for i,c in enumerate(shifted(poly,a))]
        coeffs = [sum(weights[i]*Q(comb(k,i),comb(n,i)) for i in range(k+1)) for k in range(n+1)]
        count += variations(coeffs)
    count += variations(shifted(poly,cuts[-1]))
    assert count == root_count,(count,root_count)


def cuts_for(radii, lam):
    points = [Q(0),Q(1,4),1-70*lam,Q(1),1+70*lam,Q(2)]
    for radius in radii:
        points.extend(radius*x for x in map(Q,(".1",".7",".8",".9",".99","1","1.01","2","4")))
    return sorted(set(points))


class Side:
    """A core window and optionally one remote window, in separate spatial charts."""

    def __init__(self, s, lam, even, odd, radii):
        self.profile = linear_combination([1,s],[even,odd])
        self.action = action_primitive(self.profile)
        self.windows = []
        self.center = value(even,1)
        self.lam = lam
        self.local_height = shift_at_one(self.profile,lam)
        self.local_action = shift_at_one(self.action,lam)
        local_slope = shift_at_one(deriv(self.profile),lam)
        a = -2+s
        residual = linear_combination([1/(3*lam**2),-1],[local_slope,[a,0,3]])
        for poly in (residual,deriv(residual)):
            assert sum(abs(c)*70**i for i,c in enumerate(poly)) < Q(1,100)
        k = Q(1) if a == -3 else Q(".57735")
        self.core_window = self.add_window(
            self.local_height,self.local_action,local_slope,1,1,
            (-k-Q(".02"),-k+Q(".02")),(k-Q(".02"),k+Q(".02")),
            (-2*k-Q(".03"),-2*k+Q(".03")))
        # Each successive insertion puts a well on the opposite side.
        for i,(radius,coefficient) in enumerate(radii):
            if s != (-1)**i:
                continue
            m = 8+2*i
            height_unit = coefficient*radius**m
            action_unit = height_unit**2/radius**2
            profile = scale(self.profile,radius,height_unit)
            action = scale(self.action,radius,action_unit)
            p = Q(m,m+2)
            self.add_window(profile,action,deriv(profile),height_unit,action_unit,
                            (p-Q(".01"),p+Q(".01")),(Q(".99"),Q("1.01")),
                            (Q(".001"),p-Q(".01")))
        cuts = cuts_for([r for r,c in radii],lam)
        shape_certificate(deriv(self.profile)[1:],cuts,2*len(self.windows))
        for previous,current in zip(self.windows,self.windows[1:]):
            assert previous[1][1] < current[0][0]
        # Exact global quadratic lower bound for the physical shooting test.
        lower = list(self.profile[2:])
        lower[0] -= lam**2/2
        shape_certificate(lower,cuts,0)

    def add_window(self, profile, action, slope, height_unit, action_unit,
                   peak_bracket, bottom_bracket, landing_bracket):
        peak_root = root(lambda x:value(slope,x),*peak_bracket)
        bottom_root = root(lambda x:value(slope,x),*bottom_bracket)
        bottom, peak = [poly_box(profile,r) for r in (bottom_root,peak_root)]
        assert 0 < bottom[0] <= bottom[1] < peak[0]
        lands = [root(lambda x:value(profile,x)-y,*landing_bracket) for y in bottom]
        land = poly_box(action,(lands[0][0],lands[1][1]))
        skip = minus(poly_box(action,bottom_root),land)
        assert skip[0] > 0
        self.windows.append((tuple(height_unit*x for x in bottom),
                             tuple(height_unit*x for x in peak),
                             tuple(action_unit*x for x in skip)))
        return bottom,peak,land,skip,peak_root,bottom_root

    def compact_actions(self, sample):
        target = self.center+3*self.lam**3*sample
        bottom,peak,land,skip,peak_root,bottom_root = self.core_window
        assert all(target < w[0][0] for w in self.windows[1:])
        for edge in (bottom,peak):
            assert target < edge[0] or target > edge[1]
        bracket = (-Q(70),peak_root[0]) if target < bottom[0] else (bottom_root[1],Q(70))
        entry = root(lambda x:value(self.local_height,x)-target,*bracket)
        incoming = poly_box(self.local_action,entry)
        if target > bottom[1]:
            incoming = minus(incoming,skip)
        upper = land if bottom[1] < target < peak[0] else incoming
        return upper,incoming

    def outside_action(self, radius, height_unit, target, bracket):
        profile = scale(self.profile,radius,height_unit)
        action_unit = height_unit**2/radius**2
        action = scale(self.action,radius,action_unit)
        entry = root(lambda x:value(profile,x)-target,*bracket,125)
        interval = tuple(value(action,r) for r in entry)
        for bottom,peak,skip in self.windows:
            y = target*height_unit
            if y > peak[1]:
                interval = minus(interval,tuple(x/action_unit for x in skip))
            else:
                assert y < bottom[0]
        return interval


def signed_actions(lam,even,sides,insertions):
    for i,(sample,band) in enumerate(zip(SAMPLES,WIDTHS)):
        p,n = [side.compact_actions(sample) for side in sides]
        inside(tuple(x/(9*lam**5) for x in minus(p[0],n[0])),*band)
        if i in (0,3,6,7):
            sign_interval(minus(p[1],n[1]),1 if i in (0,6) else -1)
    # The unperturbed height at z=2 is 33. Its perturbed even height is
    # a second regular sample in the same positive-sign region.
    for target in (Q(33),value(even,2)):
        data = [side.outside_action(Q(1),Q(1),target,(Q("1.99"),Q("2.01"))) for side in sides]
        sign_interval(minus(*data),1)
    for i,(radius,coefficient) in enumerate(insertions):
        m = 8+2*i
        unit = coefficient*radius**m
        t = Q(m,m+2)**(m//2)/(1+m//2)
        rows = [(Q(1,64)**m,(Q(1,128),Q(1,32)),1),
                (Q(1001,1000)*t**2,(Q(".1"),Q(2)),-1),
                (Q(8)**m*65,(Q(4),Q(16)),1)]
        for target,bracket,sign in rows:
            data = []
            for side in sides:
                # Above the peak an inner root cannot meet the target, so
                # the sign bracket in the middle row contains only the outer entry.
                data.append(side.outside_action(radius,unit,target,bracket))
            sign_interval(minus(*data),(-1)**(i+1)*sign)


def physical_reference(lam,even,odd,space,sides,scaled_launch):
    height, xmax = 50000,Q(500)
    eta = lam**2
    assert scaled_launch*space**2/eta > height
    assert all(scaled_launch > peak[1] for side in sides for bottom,peak,skip in side.windows)
    assert 2*(2*height+6) < xmax**2
    even_error = sum(i*abs(c)*(xmax/space)**(i-2)/eta for i,c in enumerate(even) if i>2)
    odd_error = sum(abs(c)*(xmax/space)**(i-3)/eta for i,c in enumerate(odd) if i>3)
    assert even[2] == eta and even_error < 1
    assert odd[3] == 5*eta and odd_error < 1
    small = Q(1,4)
    assert Q(3,2)*small**2+small+6*small**3/space < Q(1,2)
    assert 6*xmax**2/space < Q(1,2)
    memory = Q(2*height+6,2**(height-1))
    assert memory < 8*small**3/(27*space) and memory < small/27


def physical_negative(lam, even, odd, sides, radii, space):
    eta = lam**2
    eps = (eta/space)**2
    cuts = cuts_for(radii,lam)
    degree = len(even)-1
    lower_tail = even[-1]/10
    tail_start = 2*max((Q(2),*radii))
    for side in sides:
        residual = list(side.profile)
        residual[-1] -= lower_tail
        assert all(x > 0 for x in shifted(residual,tail_start))
    for radius,d,sign in zip((Q(2),*radii),(8,10,12),(1,-1,1)):
        w = value(even,radius)
        mag = sum(abs(c)*radius**i for i,c in enumerate(odd))
        numer = [Q(0)]+scale(odd,radius,mag)
        profiles = [scale(side.profile,radius,w) for side in sides]
        # Beyond 2R both profiles have a common lower bound in this chart.
        floor_coefficient = Q(1,10**13)
        for side in sides:
            lower = list(side.profile)
            lower[d] -= floor_coefficient*w/radius**d
            shape_certificate(lower,sorted({2*radius,*[x for x in cuts if x>2*radius]}),0)
        end = 2**26
        tail = w/(mag*floor_coefficient*Q(end)**(d-2)*(d-2))
        grid = [Q(0),Q(1,64)]
        while grid[-1] < end:
            start = grid[-1]
            count = 32 if start < 16 else 4
            grid.extend(start*(1+Q(i,count)) for i in range(1,count+1))
        total = (Q(0),Q(0))
        for left,right in zip(grid,grid[1:]):
            low = high = Q(1)
            for p in profiles:
                bounds = poly_box(p,(left,right))
                low *= max(Q(1),1+bounds[0])
                high *= 1+bounds[1]
            term = times(poly_box(numer,(left,right)),((right-left)/high,(right-left)/low))
            denom = 2**100
            total = plus(total,(Q((term[0]*denom).__floor__(),denom),
                                Q((term[1]*denom).__ceil__(),denom)))
        margin = min(sign*t for t in total)-tail
        assert margin > Q(1,1000)
        # w+G >= k(1+z^degree) over the whole ray. Bound the finite-epsilon
        # and physical linear corrections to the Volterra integral separately.
        k = min(w/(1+tail_start**degree),lower_tail/2)
        relative = (eta/space)/k+2*eps/k**2
        assert relative < Q(1,2) and 4*eps/k**2 < 1
        assert 4*relative/k < 2*margin*radius**2*mag/w**2


def ledger():
    upper = [1,-1]*4+[1]
    pos,neg = [-1,1,-1,1,-1],[-1]
    for t in range(310):
        m,j = 8+2*t,7+2*t
        beta = (-1)**t
        assert (upper[-1],pos[-1],neg[-1]) == (beta,-beta,-beta)
        for extra in (0,1):
            u = upper+([-beta] if extra else [])
            p = pos+([beta] if extra else [])
            q = changes(u+[-u[-1]])
            incoming = p[::-1]+neg+([beta] if extra else [])
            r = changes(incoming)
            direct = -beta if extra else beta
            assert incoming[-1] != direct
            n = max(m+2*extra+2,j+4)
            count = q+r+2
            assert (q,r) == (9+3*t+extra,4+4*t+2*extra)
            assert (n,count) == (11+2*t+extra,15+7*t+3*extra)
            assert 2*count >= 7*n-48
            if n >= 14:
                h, e = divmod(n - 10, 4)
                remote = 14*h + (9, 10, 13, 14)[e]
                assert count - remote == (2, 5, 5, 8)[e]
        signs = [-beta,beta,-beta]
        upper.extend(signs)
        pos.extend(-s for s in signs)
        neg.append(beta)
    print("The two parity counts agree through degree 630. Their gains over Proposition 5.1 are 2, 5, 5 and 8 in its four residue classes.")


def main():
    if not __debug__:
        sys.exit("Run with assertions enabled.")
    endpoint_table()
    figure_action()
    lam = Q(1,10**12)
    space = Q(10**10000)
    even,odd = balanced_core(lam)
    core_sides = [Side(s,lam,even,odd,[]) for s in (1,-1)]
    signed_actions(lam,even,core_sides,[])
    assert 100 > max(Q(33),value(even,2),value(even,1)+3*lam**3*max(SAMPLES))
    physical_reference(lam,even,odd,space,core_sides,Q(100))
    physical_negative(lam,even,odd,core_sides,[],space)
    print("The octic passes its moments, folds and actions, with a matched reference and a true negative sign.")
    insertions = []
    a = 2-Q(1,10**10)
    radii = (Q(10**90),Q(10**200))
    for i,radius in enumerate(radii):
        m,c = len(even)-1,even[-1]
        insertions.append((radius,c))
        odd = linear_combination([1,1],[odd,polynomial({m+1:(-1)**(i+1)*a*c/radius})])
        even = linear_combination([1,1],[even,polynomial({m+2:c/radius**2})])
    sides = [Side(s,lam,even,odd,insertions) for s in (1,-1)]
    signed_actions(lam,even,sides,insertions)
    print("One degree-twelve input retains all core signs and has both triples of outer-action signs.")
    far_height = 65*8**10*insertions[-1][1]*radii[-1]**10
    assert far_height > 65*8**8*insertions[0][1]*radii[0]**8 > 100
    physical_reference(lam,even,odd,space,sides,2*far_height)
    physical_negative(lam,even,odd,sides,radii,space)
    print("Three actual negative incoming signs at L=10^10000 are -,+,-, with a shootable reference.")
    ledger()


if __name__ == "__main__":
    main()
