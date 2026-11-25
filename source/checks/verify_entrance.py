"""Enclose translated and unshifted differences of incoming sextic half-orbits.

The input ordinate is the limit as x tends to positive infinity. The two
reflected fields have ordinate shift 3/500 or zero. All ODE
and sign decisions use the rational defect bounds in verify_seed.py.
"""

from fractions import Fraction as Q
import sys

from verify_seed import HIGHER, LINEAR, return_height


SHIFT = Q(3, 500)
TAIL = Q(1, 448)


def verify_tail():
    # For x >= 2 the positive even square term can be discarded. Bound each
    # possibly negative coefficient divided by x**(6-j) at the left end.
    a2, a3, a4, a5, a6 = HIGHER
    assert a2 > 0 and a4 < 0
    lower = a6 - abs(a5)/2 + a4/4 - abs(a3)/8 - max(map(abs, LINEAR))/32
    assert lower > 8
    # An inward launch below y=1 obeys 0 < dy/dx <= 1/(7*x**5).
    assert TAIL == Q(1, 7 * 4 * 2**4)


def check_ordinate(y, lower_target, upper_target, shift=SHIFT):
    bounds = []
    steps = []
    for upper in (False, True):
        # The plus-minus difference decreases with the linear coefficient.
        # Each individual terminal magnitude decreases with the ordinate on
        # x=2. Use opposite edges of the tail boxes for the two fields.
        a = LINEAR[0 if upper else 1]
        polynomial = (Q(0), a, *HIGHER)
        reflected = tuple((-1)**i * c for i, c in enumerate(polynomial))
        plus_y = y - (TAIL if upper else 0)
        minus_y = y - shift - (0 if upper else TAIL)
        plus, plus_steps = return_height(polynomial, Q(2), plus_y)
        minus, minus_steps = return_height(reflected, Q(2), minus_y)
        bounds.append(plus[1] - minus[0] if upper else plus[0] - minus[1])
        steps.extend((plus_steps, minus_steps))
    assert Q(lower_target) < bounds[0] <= bounds[1] < Q(upper_target)
    assert Q(lower_target) * Q(upper_target) > 0
    print(f"shift={shift}, y={str(y):>3}  ({lower_target}, {upper_target})   "
          f"{', '.join(map(str, steps))} accepted Taylor steps", flush=True)


if __name__ == '__main__':
    if not __debug__:
        sys.exit('Run Python with assertions enabled.')
    verify_tail()
    for row in ((Q(4, 5), '-1.1e-7', '-5e-8'),
                (Q(0), '0.003', '0.008'),
                (Q(-20), '-0.007', '-0.001')):
        check_ordinate(*row)
    # The separate exterior transfer of the decic uses these unshifted signs.
    check_ordinate(Q(1), '-2.4e-9', '-1.5e-9', shift=Q(0))
    check_ordinate(Q(0), '0.008', '0.013', shift=Q(0))
