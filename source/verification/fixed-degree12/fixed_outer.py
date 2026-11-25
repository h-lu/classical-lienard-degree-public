"""A single rational variant with an outer hyperbolic cycle."""

from fractions import Fraction as Q
import json
from pathlib import Path

from compact import coefficients, FIELD_ERROR, verify_compact
from field import field, LAMBDA, LINEAR, SEED
from outer import cap_residual, EXIT_MARGIN, hypotheses, K, peak_lift
from rational import Box, at, deriv, norm
from shooting import bound, finite_inner_geometry, seed


D = Q(6, 5)
THETA = Q('-1e-20')
A = sum(LINEAR)/2
OUTER_PATH = Path(__file__).resolve().parent/'inputs/outer-variant-field.json'


def rational_outer_interval(model):
    rho = LAMBDA**4
    reference = at(model.scaled_poly(1), D)
    return (model.height*(reference-rho),
            model.height*(reference+peak_lift(LAMBDA)+2*rho))


def verify_outer_cycle():
    model = field(D, THETA, A)
    data = json.loads(OUTER_PATH.read_text())
    assert model.physical_poly() == list(map(Q, data['coefficients']))
    h_lower, h_upper = rational_outer_interval(model)
    assert (h_lower, h_upper) == tuple(map(Q, data['outer_upper_axis_interval']))
    assert h_lower > Q('4.6') and h_lower > Q('0.012')*model.height
    assert h_upper < Q('0.014')*model.height

    enclosed = field(Box(D), Box(THETA), Box(A))
    hypotheses(enclosed, D, LAMBDA)
    finite_inner_geometry(enclosed, Box(D))
    assert cap_residual(enclosed, D, THETA, LAMBDA).strictly_inside('-7e-100', '-6.99e-100')
    # The larger minus cap provides a vertical wall for all plus section hits.
    cap_high = h_upper/model.height
    assert at(model.scaled_poly(-1), D) > cap_high
    for side in (1, -1):
        p, correction = enclosed.scaled_poly(side), enclosed.scaled_correction(side)
        assert (at(enclosed.leading_even, Q(1, 2))-norm(correction, D)).lower > Q('0.001')
        assert 6*SEED[-1]*Q('0.001')/D**6 > Q('0.01')
        assert (100*enclosed.length*correction[1]).radius_bound()+norm(deriv(correction)[1:], D) < Q('1e-8')
        assert norm(deriv(p)[1:], D) < K-Q('1e-8')

        compact_seed = seed(side, A)
        assert norm(compact_seed, Q('0.01'))+FIELD_ERROR < Q('0.002')
        assert norm(compact_seed[1:], Q('0.01'))+FIELD_ERROR < Q('0.2')
        assert norm(deriv(compact_seed), 3)+FIELD_ERROR < 10**5
        assert (enclosed.height*at(p, Q('0.01')/enclosed.length)).upper < 8
        assert Q('0.002')/model.height < Q('0.001')

    upper_at_nine = bound(1, A, 2, 8)[1]-bound(-1, A, 0, 8)[0]
    assert upper_at_nine < Q('-2e-32')
    assert bound(-1, A, 0, 8)[0] > Q('0.03')
    assert bound(-1, A, 2, 8)[1] < Q('0.1')
    # Norm estimates for the two actual divergence integrals at a crossing.
    eps, height = model.epsilon, model.height
    assert eps*(Q('0.01')/model.length)**2/(2*Q('0.009')) < Q('0.002')
    compact_time = (Q(9)+Q('0.1'))/Q('0.01') + 5
    assert compact_time < 1000
    plus_upper = (-Q('1e-5')+4*K*LAMBDA**4+2*K*eps/EXIT_MARGIN
                  +10*K/height+10**8*eps)
    minus_lower = -2*K*eps/Q('0.009')-Q('1e-8')*Q('0.014')-10**8*eps
    assert plus_upper < Q('-9e-6') and minus_lower > Q('-1e-6')
    assert abs(A) < 2 and model.physical_poly()[0] == 0 and model.physical_poly()[-1] < 0
    print('Fixed rational variant with theta=-1e-20 has |a|<2 and degree twelve.')
    print('Outer plus section launches at y=9 and y=-1 give actual D signs - and +.')
    print('At every zero between them eps*(I_+ - I_-) < -8e-6, giving a repelling cycle.')
    print('The rational upper-axis enclosure and exact coefficients are in inputs/outer-variant-field.json.')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run Python with assertions enabled.')
    fixed_fields = coefficients(OUTER_PATH)
    verify_outer_cycle()
    verify_compact(fixed_fields)
    print('Six disjoint hyperbolic limit cycles certified for this one rational variant.')
