"""Finite outer passage certificates for the degree-twelve family.

The differential inequalities that connect these tests to actual hits
and gaps are in supplement Appendix B and PROOF.md. All verification decisions use integer or
rational arithmetic. Decimal conversion is for display only.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as Q
import json
from pathlib import Path

from field import field, LAMBDA, LINEAR, SHIFT
from rational import (Box, add, at, decreasing_zero, deriv, mul, norm,
                      positive_on, primitive, scale, sub, translate)


K = 10**7
LEFT_CAP, RIGHT_CAP = Q('1.18'), Q('1.22')
EXIT, EXIT_WALL = Q('1.10'), Q('1.11')
EXIT_MARGIN = Q('1e-4')
FALLING_MARGIN = Q('0.01')
CAP_CURVATURE = Q('0.1')
TAIL = Q('0.001')


def display(q):
    with localcontext() as context:
        context.prec = 7
        return str(Decimal(q.numerator) / Decimal(q.denominator))


def hypotheses(model, d, lam):
    """The uniform polynomial and scale hypotheses for the passage proof."""
    eta, rho = lam**2, lam**4
    eps = Box.wrap(model.epsilon)
    assert Q('1e-800') < eps.lower and eps.upper < rho**4 / 8
    assert 6000*eps.upper < K*rho
    assert 4*K*rho < Q('5e-6')
    assert 20*rho < Q('0.01') and 200*K*rho < 1
    assert 2*eps.upper/Q('5e-6') < (EXIT_WALL**2-EXIT**2)/2
    assert 2*eps.upper/rho < rho
    assert rho < EXIT_MARGIN
    assert 1 < Box.wrap(model.beta).lower < Box.wrap(model.beta).upper < 5
    assert 1 < Box.wrap(model.length).lower
    assert 100 / Box.wrap(model.height).lower < rho
    assert (2*Box.wrap(model.height)/Box.wrap(model.length)**2).lower > 1

    for side in (1, -1):
        p = model.scaled_poly(side)
        perturbation = model.scaled_correction(side)
        dp = deriv(p)
        curvature = sub([0] + deriv(dp), dp)  # z**3 * phi''(u)
        assert Box.wrap(p[0]).radius_bound() == 0
        assert Box.wrap(p[12]).upper < 0
        assert norm(p, 2) + 2 < K
        assert 4*norm(dp, 2) < K
        assert 8*norm(curvature, 2) < K
        assert at(dp, LEFT_CAP).lower > 0 and at(dp, RIGHT_CAP).upper < 0
        positive_on(sub(scale(curvature, -1), [0, 0, 0, CAP_CURVATURE]), LEFT_CAP, RIGHT_CAP)
        positive_on(sub(scale(dp, -1), [0, FALLING_MARGIN]), RIGHT_CAP, Q(2))
        assert at(p, RIGHT_CAP).lower > Q('0.01')
        assert at(p, Q(2)).upper < -2
        assert all(c.upper < 0 for c in translate(add(p, [0]*12 + [TAIL]), Q(2)))
        perturbation_height = norm(perturbation, 2)
        assert at(model.leading_even, RIGHT_CAP).lower > Q('0.01')
        assert at(model.leading_even, Q(1, 2)).lower-perturbation_height > Q('0.001')
        assert perturbation_height < Q('0.001')
        # E increases up to d. Compare at d without losing its exact derivative zero.
        at_d = at(p, d)
        assert at_d.strictly_inside('0.012', '0.013')
        for z, margin in ((LEFT_CAP, Q(0)), (EXIT_WALL, 2*EXIT_MARGIN)):
            assert at_d.lower > at(model.leading_even, z).upper+perturbation_height+margin
        assert (at(deriv(perturbation), d)/d).radius_bound() < 10**5*eta
    assert Q('0.013') + peak_lift(lam) < Q('0.014')
    assert LEFT_CAP < Box.wrap(d).lower <= Box.wrap(d).upper < RIGHT_CAP


def peak_lift(lam):
    return (10**5*lam**2)**2 / (2*CAP_CURVATURE)


def cap_residual(model, d, theta, lam):
    """Enclose the difference of actual backwards hits from t=-1."""
    odd_at_d = lam**2*model.beta*theta + lam**3*at(model.cubic_odd, d) + at(model.odd_seed, d)
    error = peak_lift(lam) + 3*lam**4
    return 2*odd_at_d + Box(-error, error)


def falling_action(model, side, ordinate):
    p = model.scaled_poly(side)
    dp = deriv(p)
    peak = decreasing_zero(lambda z: at(dp, z)/z, LEFT_CAP, RIGHT_CAP)
    root = decreasing_zero(lambda z: at(p, z)-ordinate, RIGHT_CAP, Q(2))
    squared = mul(dp, dp)
    integrate = primitive(squared[1:])
    # The omitted term is p[1]**2 * log(root / peak), between 0 and p[1]**2.
    return at(integrate, root) - at(integrate, peak) + Box(0, squared[0].upper)


def gap_cost(model, side, ordinate, lam):
    """Enclose -eps*log(H_s^* - h_s), for the launches in PROOF.md."""
    error = 10*K*lam**4
    return falling_action(model, side, ordinate) + Box(-error, error)


def terminal_checks(model, d, lam):
    """Uniform final lower-axis sample and the size of the matching gaps."""
    eta = lam**2
    eps = model.epsilon
    section_cost = gap_cost(model, 1, Box(-22, 1)/model.height, lam)
    direct_costs = [gap_cost(model, s, -Q(1, 2), lam) for s in (1, -1)]
    assert section_cost.strictly_inside('0.0019', '0.0021')
    for cost in direct_costs:
        assert cost.strictly_inside('1.3', '1.4')
        assert (cost/eps).lower > Q('3e780')
    # Both section-gap costs obey the same bound. Certify the second side
    # separately for the necessity of cap agreement at physical sections.
    other_cost = gap_cost(model, -1, Box(-22, 1)/model.height, lam)
    for cost in (section_cost, other_cost):
        assert cost.strictly_inside('0.0019', '0.0021')
        assert (cost/eps).lower > Q('4e777')

    # Compare the direct action difference without subtracting interval
    # copies of the shared even action over the d-box.
    even = model.leading_even
    even_slope = deriv(even)[1:]
    odd = scale(model.odd, 1/eta)
    wronskian = sub(mul(even_slope, deriv(odd)), mul(deriv(even_slope), odd))
    root = decreasing_zero(lambda z: at(even, z)+Q(1, 2), RIGHT_CAP, Q(2), steps=100)
    integrate = primitive(wronskian)
    first_order = at(integrate, root)-at(integrate, d)
    perturbation_bound = 10**7*eta
    # On the homotopy from the even profile, each second derivative of the
    # action in the homotopy parameter is bounded by this remainder.
    for side in (1, -1):
        p = model.scaled_correction(side)
        assert max(norm(p, 2), norm(deriv(p), 2)) < perturbation_bound
    assert at(deriv(even), LEFT_CAP).lower > 0
    assert at(deriv(even), RIGHT_CAP).upper < 0
    assert at(even, RIGHT_CAP).lower > Q('0.01') and at(even, Q(2)).upper < -2
    assert 4*norm(deriv(even), 2) < K
    assert 8*norm(sub([0]+deriv(deriv(even)), deriv(even)), 2) < K
    positive_on(sub(scale(deriv(even), -1), [0, FALLING_MARGIN]), RIGHT_CAP, Q(2))
    positive_on(sub(scale(sub([0]+deriv(deriv(even)), deriv(even)), -1),
                    [0, 0, 0, CAP_CURVATURE]), LEFT_CAP, RIGHT_CAP)
    error = 10**5*K*perturbation_bound**2 + 20*K*lam**4
    residual = -2*first_order + Box(-error, error)/eta
    assert residual.strictly_inside('3.8', '5.5')
    print('  At a joint match the direct t=-1/2 sample has actual return sign +.')
    print('  Uniform (eps/eta)*log(g_+^b/g_-^b) in (3.8, 5.5).')
    print('  On the box: section gaps < exp(-4e777), direct gaps < exp(-3e780).')


def second_scale():
    """One prescribed larger splitting scale, with no parameter search."""
    lam, d = Q('1e-10'), Q(6, 5)
    model = field(Box(d), Box(0), Box(sum(LINEAR)/2), lam)
    hypotheses(model, d, lam)
    for side in (1, -1):
        cost = gap_cost(model, side, 0, lam)
        assert (cost/model.epsilon).strictly_inside('4.81e177', '4.83e177')
    print('At lambda=1e-10, d=6/5, theta=0, a=mid J: true -log(g_s(0)) in (4.81e177, 4.83e177).')


def check_diagnostic():
    path = Path(__file__).resolve().parent / 'inputs/fixed-diagnostic-field.json'
    data = json.loads(path.read_text())
    d, theta, lam = (Q(data[key]) for key in ('d', 'theta', 'lambda'))
    a = sum(LINEAR)/2
    model = field(d, theta, a, lam)
    assert model.physical_poly() == list(map(Q, data['coefficients']))
    assert tuple(map(Q, (data['L'], data['S'], data['epsilon']))) == (model.length, model.height, model.epsilon)
    enclosed = field(Box(d), Box(theta), Box(a), lam)
    hypotheses(enclosed, d, lam)
    delta_h = cap_residual(enclosed, d, theta, lam)
    assert delta_h.strictly_inside('6.99e-100', '7e-100')
    # The rational midpoint of the two heights P_s(d) lies between the
    # smaller true peak (plus 2*rho) and the larger actual reference hit.
    radius_midpoint = at(add(model.leading_even, model.even_correction), d)
    odd_at_d = at(model.odd, d)
    assert odd_at_d > peak_lift(lam)+2*lam**4
    assert radius_midpoint > 0
    assert (model.height*radius_midpoint) > Q('8e465')
    assert 2*model.epsilon/lam**4 < lam**4
    assert 2*model.epsilon < Q('0.001')
    print('Fixed diagnostic polynomial equals its exact rational input.')
    print('Actual H_+^* - H_-^* in (6.99e-100, 7e-100), scaled upper-axis heights.')
    print('At h=S*(Phi_+(d)+Phi_-(d))/2 the plus half returns and the minus half escapes.')
    print('Rational escape-test h approximately ' + display(model.height*radius_midpoint))

    gaps = []
    for side in (1, -1):
        gap = gap_cost(enclosed, side, Box(-20, 1)/enclosed.height, lam)
        assert gap.strictly_inside('0.0019848', '0.0019850')
        exponent = gap / enclosed.epsilon
        assert exponent.strictly_inside('4.81e777', '4.83e777')
        gaps.append(gap)
    log_ratio = (gaps[1]-gaps[0]) / lam**2
    assert log_ratio.strictly_inside('8.57e-5', '8.58e-5')
    print('Actual -log(g_s(y/S)) in (4.81e777, 4.83e777) for both sides and -20 <= y <= 1.')
    print('Actual (eps/eta)*log(g_+(Delta/S)/g_-(0)) in (8.57e-5, 8.58e-5).')


def check_matching_box(lam, radii, theta_radius):
    linear = Box(*LINEAR)
    d, theta = Box(*radii), Box(-theta_radius, theta_radius)
    model = field(d, theta, linear, lam)
    hypotheses(model, d, lam)
    print(f'Outer matching box at lambda={display(lam)}, theta radius={display(theta_radius)}')
    for sign, interval in ((-1, ('-46.75', '-46.74')), (1, ('92.96', '92.97'))):
        theta_face = sign*theta_radius
        face = field(d, Box(theta_face), linear, lam)
        residual = cap_residual(face, d, theta_face, lam) / lam**3
        assert residual.strictly_inside(*interval)
        print(f'  Theta face {sign:+d} gives (H_+^*-H_-^*)/lambda^3 in {interval}.')
    for sign, endpoint, interval in zip((-1, 1), radii,
                                        (('-2.974e-13', '-2.972e-13'), ('4.704e-13', '4.707e-13'))):
        face = field(Box(endpoint), theta, linear, lam)
        plus = gap_cost(face, 1, SHIFT/face.height, lam)
        minus = gap_cost(face, -1, 0, lam)
        residual = (minus-plus) / lam**2
        assert residual.strictly_inside(*interval)
        print(f'  Radius face {sign:+d} gives (eps/eta)*log(g_+(Delta/S)/g_-(0)) in {interval}.')
    estimate = -at(model.cubic_odd, d)/model.beta
    error = ((2*at(model.odd_seed, d)).radius_bound()+peak_lift(lam)+3*lam**4)
    interval = estimate + Box(-error, error) / (2*lam**3*model.beta)
    assert interval.strictly_inside('-3.4', '-3.2')
    print('  At any actual cap match theta/lambda belongs to (-3.4, -3.2).')
    print(f'  At a cap match |theta + lambda*H(d)/beta| < {display(10**10*lam**2)}')
    assert (Box(error)/(2*lam**2*model.beta)).upper < 10**10*lam**2
    return model


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run Python with assertions enabled.')
    check_diagnostic()
    second_scale()
