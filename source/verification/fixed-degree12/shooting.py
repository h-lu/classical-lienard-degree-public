"""Finite uniform compact shooting faces and section entrance inequalities."""

from fractions import Fraction as Q

from compact import half_return
from field import field, LAMBDA, LINEAR, RADIUS, SEED, SHIFT
from outer import check_matching_box, K, terminal_checks
from rational import Box, add, at, deriv, norm, positive_on, sub


PERTURBATION = Q('1e-72')
ENTRANCE_SHIFT_ERROR = Q('1e-12')
INWARD_TAIL = Q(1, 448)
CAPTURE_WALL = Q('1.18')


def seed(side, a):
    return tuple(side**i*c for i,c in enumerate((Q(0), a, *SEED)))


def bound(side, a, x, y, field_error=PERTURBATION):
    return half_return(seed(side, a), Q(x), Q(y), field_error=field_error)[0]


def finite_inner_geometry(model, d):
    length, height = model.length, model.height
    assert length.lower > 20
    reference = (at(model.leading_even, 1)+at(model.leading_even, d))/2
    assert (height * reference).lower > 9
    assert (at(model.leading_even, Q(1, 2)) - reference).upper < 0
    assert (at(model.leading_even, Q('1.1')) - reference).upper < 0
    for side in (1, -1):
        poly = model.scaled_poly(side)
        assert (at(poly, CAPTURE_WALL)-reference).lower > 0

        # Coefficients below degree six agree with the seed, including the
        # exact quintic cancellation in the definition of L.
        perturbation = [Q(0)]*6 + [6*SEED[-1]*model.scaled_correction(side)[6]]
        perturbation += [6*SEED[-1]*poly[i]/length**(i-6) for i in range(7, len(poly))]
        for order in (0, 1):
            assert norm(perturbation, 3) < PERTURBATION
            perturbation = deriv(perturbation)

        # The following uniform inequalities implement the inward tail
        # comparisons at [2,10], [10,L/2] and [L/2,L*CAPTURE_WALL].
        a2, a3, a4, a5, a6 = SEED
        amax = max(abs(a) for a in LINEAR)
        assert a2 > 0 and a4 < 0
        seed_tail = a6-abs(a5)/2+a4/4-abs(a3)/8-amax/32
        assert seed_tail - Q('1e-68') > 8
        physical = [model.scaled_correction(side)[6]]
        physical += [poly[i]/length**(i-6) for i in range(7, len(poly))]
        assert 6*a6*norm(physical, 10) < Q('1e-68')
        # For x >= 10, z <= 1/2, E(z) >= z**6 (9/16)(1-1/(4*d**2))/6.
        even_lower = a6 * Q(9, 16) * (1 - 1/(4*d**2))
        remainder = norm(model.scaled_correction(side)[6:], Q(1, 2))
        assert (even_lower+a4/100-abs(a5)/10-abs(a3)/1000-amax/10**5-6*a6*remainder).lower > 8
        assert (6*a6*Q('0.001')/CAPTURE_WALL**6) > Q('0.01')
        assert (at(model.leading_even, Q(1, 2))-norm(model.scaled_correction(side), CAPTURE_WALL)).lower > Q('0.001')
    # Capture of v=x**2/2 from arbitrary v above 50, then below 2.
    assert 9 < 4*50**3/100 and 9 < (64-4)*2**3
    assert Q(25, 2*50**2) + Q(1, 32) < 1


def reference_faces():
    low, high = LINEAR
    lower = bound(1, low, 0, 8)[0]-bound(-1, low, 2, 8)[1]
    upper = bound(1, high, 2, 8)[1]-bound(-1, high, 0, 8)[0]
    assert lower > Q('2e-32') and upper < Q('-6e-32')
    print('Actual D(ST) > 2e-32 at a_- and < -6e-32 at a_+ uniformly in d, theta.')


def diagnostic_reference():
    d = Q(6, 5)
    mid = sum(LINEAR)/2
    model = field(Box(d), Box('1e-20'), Box(mid))
    finite_inner_geometry(model, Box(d))
    upper = bound(1, mid, 2, 8)[1]-bound(-1, mid, 0, 8)[0]
    assert upper < Q('-2e-32')
    print('In the original fixed diagnostic field, actual D(ST) < -2e-32.')


def compact_signs():
    heights = tuple(map(Q, ('.5', '3.7', '3.88', '4.03', '4.36')))
    margins = tuple(map(Q, ('2e-7', '1e-24', '3e-27', '2e-28', '6e-32')))
    for index, (height, margin) in enumerate(zip(heights, margins)):
        orientation = (-1)**index
        a = LINEAR[1 if orientation > 0 else 0]
        p, m = (bound(s, a, 0, height) for s in (1, -1))
        lower = p[0]-m[1] if orientation > 0 else m[0]-p[1]
        assert lower > margin
    print('Uniform actual compact signs in the joint parameter box: +,-,+,-,+.')


def shoulder_obstruction(model, d):
    """Uniform contraction before the compact section, independent of fold actions."""
    left, right = Q(1, 2), Q(9, 10)
    e = model.leading_even
    sample_base = at(add(e, model.even_correction), 1)
    unit = Q(4, 3)*(1-1/d**2)
    reference = (at(e, 1)+at(e, d))/2
    height_window = sample_base + unit*LAMBDA**3*Box(-100000, 100000)
    assert height_window.strictly_inside('0.0108', '0.0110')
    assert height_window.upper < reference.lower
    inside_error = Q('1e-8')
    for side in (1, -1):
        p, delta = model.scaled_poly(side), model.scaled_correction(side)
        perturbation_bound = norm(delta, CAPTURE_WALL)
        assert (at(e, left)+perturbation_bound).upper < Q('0.0019')
        assert (at(e, right)-perturbation_bound).lower > Q('0.0071')
        positive_on(sub(deriv(p), [0, Q('0.005')]), left, right)
        linear_term = 10*Box.wrap(delta[1]).radius_bound()*model.length.upper
        assert norm(deriv(delta)[1:], 1)+linear_term < inside_error
        for a in LINEAR:
            polynomial = seed(side, a)
            assert norm(polynomial, Q('0.1'))+PERTURBATION < 8
            assert norm(polynomial, 3)+PERTURBATION+8 < 10**5
            assert norm(deriv(polynomial), 3)+PERTURBATION < 10**5
        # Parameter monotonicity gives the required w bounds from the two faces.
        # The enlarged defect verifies the time and lower magnitude of the
        # q=0 arc uniformly in a for this terminal conversion.
        uniform_error = Q('1e-29')
        assert 3*(LINEAR[1]-LINEAR[0])+PERTURBATION < uniform_error
        assert bound(side, LINEAR[1 if side == 1 else 0], 0, 8, uniform_error)[0] > Q('0.03')
        assert bound(side, LINEAR[0 if side == 1 else 1], 2, 8)[1] < Q('0.1')
    assert Q(2)*model.epsilon.upper/Q('1e-4') < Q('0.001')
    assert Q(8)/model.height.lower < Q('0.002')
    assert Q('0.004')*Q('0.005')-Q('0.002')*inside_error > Q('1.9e-5')
    # Terminal conversion 400 L^2 < exp(1e-6/eps), by exp(t)>t.
    assert 400*model.length.upper**2*model.epsilon.upper < Q('1e-6')
    assert Q('1.8e-5')/model.epsilon.upper > Q('4e775')
    print('  Any mixed shoulder signs require |D(ST)| < exp(-4e775).')
    print('  For fixed d,theta this requires |a-a_*(d,theta)| < exp(-3e775).')


def fine_section_checks(model):
    eps = model.epsilon.upper
    height_lo, height_hi = model.height.lower, model.height.upper
    assert Q('1e-800') < model.epsilon.lower and 2**3000 > 10**800
    assert 25/height_lo + 400000*eps < Q('1e-4')
    assert 1500*eps < Q('1e-4')
    assert 4*K*LAMBDA**4 < Q('1e-6')
    assert 12*eps/Q('2e-6')**3 < 1
    assert 2*(2 + 2*K*400000 + K*10002) + 5*K*eps < 3*10**13
    assert 3*10**13*eps*height_hi < Q('1e-280')
    for side in (1, -1):
        p = model.scaled_correction(side)
        # Compare both falling slopes with the leading even slope at t=0.
        slope_change = norm(deriv(p), 2)+100*K*(norm(p, 2)+25/height_lo)
        assert 26*slope_change < Q('1e-55')/2
        assert (at(model.leading_even, Q('1.22'))).lower > Q('0.01')
    assert Q('0.01')*ENTRANCE_SHIFT_ERROR > 2*Q('1e-55')
    print('At matching hits, fine entrance pairing obeys |y_minus - (y_plus-3/500)| < 1e-12.')


def section_signs():
    for ordinate, orientation, margin in ((Q(4, 5), -1, Q('5e-8')),
                                           (Q(0), 1, Q('0.003')),
                                           (Q(-20), -1, Q('0.001'))):
        a = LINEAR[1 if orientation > 0 else 0]
        if orientation > 0:
            wplus = bound(1, a, 2, ordinate)
            wminus = bound(-1, a, 2, ordinate-SHIFT-ENTRANCE_SHIFT_ERROR-INWARD_TAIL)
            lower = wplus[0]-wminus[1]
        else:
            wplus = bound(1, a, 2, ordinate-INWARD_TAIL)
            wminus = bound(-1, a, 2, ordinate-SHIFT+ENTRANCE_SHIFT_ERROR)
            lower = wminus[0]-wplus[1]
        assert lower > margin
    print('At any joint match, actual outer section return signs: -,+,-.')


def main():
    model = check_matching_box(LAMBDA, RADIUS, 10*LAMBDA)
    assert model.epsilon.strictly_inside('4.07e-781', '4.09e-781')
    finite_inner_geometry(model, Box(*RADIUS))
    reference_faces()
    diagnostic_reference()
    compact_signs()
    shoulder_obstruction(model, Box(*RADIUS))
    fine_section_checks(model)
    section_signs()
    terminal_checks(model, Box(*RADIUS), LAMBDA)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run Python with assertions enabled.')
    main()
