"""Exact arithmetic for the primary-source comparisons.

DMH, Proposition 2 and §3.1, doi:10.1007/s10884-014-9358-1.
Odani, Example 3, corrected factorization on p. 230, doi:10.24546/0100499900.
Palomar, Appendix A, paper/main.tex at commit 7da225327c6c.
OpenAI quintic manuscript, §7, at commit adc7f1241b42.
"""

from fractions import Fraction as Q
import sys

from verify_outer import value
from verify_seed import HIGHER, LINEAR, multiply
from verify_square_shoulders import deriv, linear_combination, primitive


# Appendix A, rational sextic center, in order a_1 through a_6.
PALOMAR_COEFFICIENTS = (
    Q("-92025883280136345928114/26103339339473830880164235"),
    Q("259235225707582224718121441/18715972094399964530982824"),
    Q("10672739305352820569952968/12336204592094761858066741"),
    Q("-336421133612586705537064207/12229582503151462140066735"),
    Q("-12665052643337593541339158/7109329057417444170418853"),
    Q("233599035036509246393129479/11502661198404994310442156"),
)
# Appendix A rounds each terminal magnitude outward to 42 decimal places.
# Each row gives lower/upper bounds for the positive side, then the reflected side.
PALOMAR_ENDPOINTS = (
    ("0.036491285574418155644860079650087057000070", "0.036491285574418155644860079650087057046785",
     "0.036491069951363189750046519718097492326113", "0.036491069951363189750046519718097492372804"),
    ("0.036491804793302773054131196805446637178353", "0.036491804793302773054131196805446637293390",
     "0.036491804793302773054134115535427368812008", "0.036491804793302773054134115535427368925955"),
    ("0.036491804793302773060509942667524930455231", "0.036491804793302773060509942667524930572526",
     "0.036491804793302773060509939340402927841001", "0.036491804793302773060509939340402927957760"),
    ("0.036491804793302773060737445440385434619820", "0.036491804793302773060737445440385434738972",
     "0.036491804793302773060737445689481225676454", "0.036491804793302773060737445689481225795425"),
    ("0.036491804793302773060752290676915195030000", "0.036491804793302773060752290676915195153478",
     "0.036491804793302773060752290676778157621737", "0.036491804793302773060752290676778157745014"),
    ("0.036491804793302773060752304184606198227803", "0.036491804793302773060752304184606198353873",
     "0.036491804793302773060752304184608498616025", "0.036491804793302773060752304184608498742093"),
    ("0.036491804793302773060752304327613676693144", "0.036491804793302773060752304327613676820892",
     "0.036491804793302773060752304327611376304952", "0.036491804793302773060752304327611376432642"),
)
ABSOLUTE_GAPS = ("2.15e-7", "2.91e-24", "3.32e-27", "2.49e-28",
                 "1.37e-31", "2.30e-33", "2.30e-33")


def slow_divergence_sextic():
    alpha = list(map(Q, (1, 0, "-3.1", 0, "2.7")))
    beta = list(map(Q, (0, 0, 0, "-0.4", 0, 1)))
    assert alpha[2]**2 < 4*alpha[0]*alpha[4]
    integrated = primitive(linear_combination(
        (1, -1), (multiply(deriv(alpha), beta), multiply(alpha, deriv(beta)))))
    cubic = [Q(2, 5), -Q(156, 125), Q(411, 350), -Q(3, 10)]
    assert integrated == [0, 0, 0, cubic[0], 0, cubic[1], 0, cubic[2], 0, cubic[3]]
    values = [value(cubic, u) for u in (0, Q(3, 4), 1, 3)]
    assert values == [Q(2, 5), -Q(227, 112000), Q(23, 875), -Q(766, 875)]
    print("De Maesschalck-Huzak sextic: three simple positive cubic roots by rational signs.")


def palomar_intervals():
    assert len(PALOMAR_ENDPOINTS) == len(ABSOLUTE_GAPS) == 7
    for i, (row, gap) in enumerate(zip(PALOMAR_ENDPOINTS, ABSOLUTE_GAPS)):
        plus_lo, plus_hi, minus_lo, minus_hi = map(Q, row)
        assert plus_lo < plus_hi and minus_lo < minus_hi
        lower_signed_gap = plus_lo-minus_hi if i % 2 == 0 else minus_lo-plus_hi
        assert lower_signed_gap > Q(gap)

    assert tuple(round(v, 7) for v in PALOMAR_COEFFICIENTS[1:]) == HIGHER
    difference = abs(PALOMAR_COEFFICIENTS[1] - HIGHER[0])
    assert Q("2.75e-8") < difference < Q("2.76e-8")
    assert PALOMAR_COEFFICIENTS[0] + Q("1e-40") < LINEAR[0]
    print("Palomar Appendix A: seven alternating strict interval gaps; the rounded seed lies outside the box.")


def odani_quintic():
    def squared_primitive(u, nu):
        return u * (u**2 - nu*u + 1)**2

    # Write z=x^2 and r=(nu+c)/2, so nu=r+1/r, c=r-1/r, beta_1^2=1/r.
    # After multiplication by r^5 both sides have degree at most 10 in r
    # and at most 5 in z. The full grid therefore proves the identity.
    for r in map(Q, range(2, 13)):
        nu = r + 1/r
        c = r - 1/r
        for z in map(Q, range(6)):
            difference = squared_primitive(z+c, nu) - squared_primitive(z, nu)
            factor = c * (z - 1/r)**2 / 5
            square = (5*z - (3*nu-5*c)/2)**2
            assert difference == factor * (square + 4*(nu**2-5))
            # Page 230 of the Kobe scan, last factored RHS.
            d = (3*nu - 5*c)/2
            printed = c*(z - 1/r)**2 * (5*(z-d)**2 + 16*(nu**2-5))
            normalized_with_16 = factor * (square + 16*(nu**2-5))
            if (r, z) == (2, 0):
                assert (difference, printed, normalized_with_16) == (
                    Q(3, 8), Q(15, 2), Q(3, 2))
            if (r, z) == (3, 1):
                # Both z and the shift d are nonzero in this comparison.
                assert (difference, printed, normalized_with_16) == (
                    Q(3968, 243), Q(12800, 81), Q(8192, 243))
    print("Odani Example 3: corrected factorization holds on the identity grid.")
    print("At mu=1, nu=5/2, z=0: actual = 3/8, printed = 15/2, normalized with 16 = 3/2.")


def quintic_energy():
    # First energy displacement in §7 divided by pi.
    energy = [0, 0, -4, 0, 5, 0, -1]
    assert [value(energy, s) for s in (1, 2)] == [0, 0]
    assert [value(deriv(energy), s) for s in (1, 2)] == [6, -48]
    print("Quintic preprint: the two sharpness roots have normalized derivatives 6 and -48.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run Python with assertions enabled.")
    slow_divergence_sextic()
    odani_quintic()
    palomar_intervals()
    quintic_energy()
