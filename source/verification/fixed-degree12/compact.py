"""Validated compact half-returns and divergence integrals for the fixed field."""

from fractions import Fraction as Q
import json
from math import lcm
from pathlib import Path

from field import LINEAR, SEED
from rational import deriv, sub
from reference.checks.verify_seed import (
    SCALE, at_time, ceil_div, check_derivative, multiply, residual_norm,
    taylor_path, tube,
)


ROOT_DIR = Path(__file__).resolve().parent
DOMAIN = Q(3)
FIELD_ERROR = Q('2e-74')


def coefficients(path=ROOT_DIR/'inputs/fixed-diagnostic-field.json'):
    document = json.loads(path.read_text())
    full = list(map(Q, document['coefficients']))
    seed = [Q(0), sum(LINEAR)/2, *SEED]
    difference = sub(full, seed)
    for derivative_order in (0, 1):
        assert sum(abs(coef)*DOMAIN**i for i, coef in enumerate(difference)) < FIELD_ERROR
        difference = deriv(difference)
    return [tuple(s**i*c for i, c in enumerate(seed)) for s in (1, -1)]


def integration_grid(derivative, path, inv_step):
    """Dyadic coefficients enclosing the primitive of h*p'(X(tau))."""
    denominator = lcm(*(c.denominator for c in derivative))
    numerators = [int(c * denominator) for c in derivative]
    degree = len(numerators) - 1
    composed = [numerators[-1]]
    for j in range(degree-1, -1, -1):
        composed = multiply(path, composed)
        composed[0] += numerators[j] * SCALE**(degree-j)
    unit = denominator * SCALE**(degree-1) * inv_step
    lower, upper = [0], [0]
    for i, coefficient in enumerate(composed, 1):
        lower.append(coefficient // (unit*i))
        upper.append(ceil_div(coefficient, unit*i))
    return lower, upper


def event_bounds(xs, ys, error, x_bound, speed):
    left, right = 0, SCALE
    for _ in range(200):
        mid = (left + right) // 2
        xlo, xhi = at_time(xs, mid)
        if xlo > 0:
            left = mid
        elif xhi < 0:
            right = mid
        else:
            break
    x_error = Q(max(abs(xlo), abs(xhi))+error, SCALE)
    time_error = x_error / speed
    y_error = Q(error, SCALE) + x_bound*time_error
    ylo, yhi = at_time(ys, mid)
    result = -Q(yhi, SCALE)-y_error, -Q(ylo, SCALE)+y_error
    assert result[0] > 0
    return result, mid, time_error


def half_return(seed, x0, lower_h, upper_h=None, integral=False, field_error=FIELD_ERROR):
    """Enclose the hits for the full fixed field, using a sextic path proposal.

    An initial height interval encloses every trajectory in the same tube.
    The perturbation between the fixed field and the seed enters the defect
    and the divergence integral. Initial x0 is fixed.
    """
    if upper_h is None:
        upper_h = lower_h
    check_derivative(seed)
    denominator = lcm(*(c.denominator for c in seed))
    numerators = [int(c * denominator) for c in seed]
    derivative = deriv(seed)
    second_bound = sum(abs(c)*DOMAIN**i for i, c in enumerate(deriv(derivative)))
    x, y = (x0*SCALE).__floor__(), ((lower_h+upper_h)*SCALE/2).__floor__()
    error = (((upper_h-lower_h)*SCALE/2).__ceil__())+2
    elapsed, steps, inv_step = Q(0), 0, 16
    log_lo, log_hi = Q(0), Q(0)
    while elapsed < 40:
        assert inv_step < 2**32
        xs, ys = taylor_path(seed, inv_step, x, y)
        num, den = residual_norm(numerators, denominator, inv_step, xs, ys)
        if num * inv_step * 10**53 >= den:
            inv_step *= 2
            continue
        local = (Q(SCALE*num, den) + SCALE*field_error/inv_step).__ceil__()
        end_error = ceil_div((error+local)*inv_step, inv_step-1)
        xmin, xmax = tube(xs, end_error)
        if 100*xmin <= -SCALE or xmax >= SCALE*DOMAIN:
            inv_step *= 2
            continue
        if xmin > 0:
            # The same SOS bound proves p' >= p'(0) > -1/250 for x >= 0.
            assert seed[1] > -Q(1, 250)
            factor = 250*inv_step
            end_error = ceil_div((error+local)*factor, factor-1)
            xmin, xmax = tube(xs, end_error)
        xend, yend = sum(xs), sum(ys)
        terminal = False
        if xmin <= 0:
            ymin, ymax = tube(ys, end_error)
            x_bound = Q(max(abs(xmin), abs(xmax)), SCALE)
            flux_bound = sum(abs(c)*x_bound**j for j, c in enumerate(seed))+field_error
            departing = steps == 0 and x0 == 0 and Q(ymin, SCALE) > flux_bound
            speed = -Q(ymax, SCALE) - flux_bound
            if not departing:
                if speed <= 0 or x <= error:
                    inv_step *= 2
                    continue
                terminal = xend < -end_error
                if not terminal and xend <= end_error:
                    inv_step *= 2
                    continue
        if integral:
            floored, ceiled = integration_grid(derivative, xs, inv_step)
            local_integral_error = (second_bound*Q(end_error, SCALE)+field_error)/inv_step
        if terminal:
            result, location, time_error = event_bounds(xs, ys, end_error, x_bound, speed)
            if integral:
                max_deriv = sum(abs(c)*x_bound**j for j, c in enumerate(derivative))+field_error
                width = time_error*max_deriv+local_integral_error
                log_lo += Q(at_time(floored, location)[0], SCALE)-width
                log_hi += Q(at_time(ceiled, location)[1], SCALE)+width
            return result, (log_lo, log_hi), steps+1
        if integral:
            log_lo += Q(sum(floored), SCALE)-local_integral_error
            log_hi += Q(sum(ceiled), SCALE)+local_integral_error
        x, y, error = xend, yend, end_error
        elapsed += Q(1, inv_step)
        steps += 1
        if inv_step > 8 and num * inv_step * 10**66 < den:
            inv_step //= 2
    raise ArithmeticError('return time exceeded 40')


def return_difference(fields, height):
    p, m = (half_return(seed, Q(0), Q(height))[0] for seed in fields)
    return p[0]-m[1], p[1]-m[0]


def divergence_difference(fields, lower, upper):
    p, m = (half_return(seed, Q(0), lower, upper, integral=True)[1] for seed in fields)
    return p[0]-m[1], p[1]-m[0]


def sign(box):
    if box[0] > 0:
        return 1
    if box[1] < 0:
        return -1
    raise ArithmeticError('uncertain sign')


def verify_compact(fields):
    path = ROOT_DIR/'inputs/compact-intervals.json'
    intervals = json.loads(path.read_text())
    assert len(intervals) == 5
    integral_bounds = (('0.00235', '0.00236'), ('-4.1e-5', '-3.6e-5'),
                       ('1.5e-5', '2.1e-5'), ('-1.5e-5', '-9e-6'), ('1.8e-5', '2.6e-5'))
    previous = Q(0)
    print('Fixed rational degree-12 field, exact compact half-returns and divergence.')
    for index, (left, right) in enumerate(intervals):
        lo, hi = Q(left), Q(right)
        assert previous < lo < hi
        previous = hi
        expected = (-1)**index
        assert sign(return_difference(fields, lo)) == expected
        assert sign(return_difference(fields, hi)) == -expected
        log_difference = divergence_difference(fields, lo, hi)
        assert sign(log_difference) == expected
        lower, upper = integral_bounds[index]
        assert Q(lower) < log_difference[0] <= log_difference[1] < Q(upper)
        print(f'  ({left}, {right}): signs {expected:+}/{-expected:+}, '
              f'integral difference in ({lower}, {upper})',
              flush=True)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run Python with assertions enabled.')
    verify_compact(coefficients())
