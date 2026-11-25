"""Degree-twelve family from Section S5 of the electronic supplement."""

from dataclasses import dataclass
from fractions import Fraction as Q

from rational import add, at, deriv, mul, primitive, scale, sub


SEED = tuple(map(Q, ('13.8510158', '0.8651558', '-27.5087995', '-1.7814695', '20.3082601')))
LINEAR = tuple(map(Q, ('-0.003525444747007079096354237286',
                       '-0.003525444747007079096354237285')))
RADIUS = tuple(map(Q, ('1.19888999543', '1.19888999544')))
LAMBDA = Q('1e-40')
SHIFT = Q(3, 500)


@dataclass
class Field:
    leading_even: list
    even_correction: list
    odd: list
    beta: object
    cubic_odd: list
    odd_seed: list
    length: object
    height: object
    epsilon: object

    def scaled_poly(self, side):
        return add(self.leading_even, self.scaled_correction(side))

    def scaled_correction(self, side):
        return add(self.even_correction, scale(self.odd, side))

    def physical_poly(self):
        return [self.height * c / self.length**i for i, c in enumerate(self.scaled_poly(1))]


def moment(a, p):
    return at(primitive(mul(a, deriv(p))), 1)


def field(d, theta, linear, lam=LAMBDA):
    """Use Fraction for a point and Box for interval coefficients."""
    eta = lam**2
    shoulder = [1, 0, -1]
    cap = [1, 0, -1/d**2]
    both = mul(shoulder, cap)
    a = [0]*4 + mul(shoulder, both)
    even = primitive([0] + a)
    weight = [0]*8 + mul(both, both)
    k = -at(primitive(mul(weight, [1, 0, 1])), 1) / at(primitive(mul(weight, [0, 0, 3, 0, -1])), 1)
    base = scale([0]*5 + mul(both, [1, 0, k]), -1)
    u = [0]*7 + shoulder
    u_cap = mul(u, cap)
    cap_mode = scale(sub(u, scale(u_cap, moment(a, u) / moment(a, u_cap))), 1/at(u, d))
    unscaled = add(base, scale(cap_mode, theta))
    unit = Q(4, 3) * (1 - 1/d**2)
    beta = -unit / at(deriv(unscaled), 1)
    odd_linear = scale(unscaled, beta)
    odd_cubic = [0]*7 + [8*unit]
    odd_cubic = sub(odd_cubic, scale(u_cap, moment(a, odd_cubic)/moment(a, u_cap)))
    even_mode = odd_linear[1:]
    pivot = moment(even_mode, odd_linear)
    control = [0]*4 + [-2*unit]
    control = sub(control, scale(even_mode, moment(control, odd_linear)/pivot))
    corrected_control = scale(even_mode, (unit**2*Q(493, 8) - moment(control, odd_cubic))/pivot)
    length = -SEED[-2]/(6*SEED[-1]*beta*eta)
    height = 6*SEED[-1]*length**6
    epsilon = length**2 / height**2
    seed = [0, *(coefficient * length**i / height
                 for i, coefficient in enumerate((linear, *SEED[:3]), 1))]
    odd_seed = [0, seed[1], 0, seed[3]]
    even_correction = add(primitive([0] + add(scale(control, eta), scale(corrected_control, lam**3))),
                          [0, 0, seed[2], 0, seed[4]])
    odd = add(add(scale(odd_linear, eta), scale(odd_cubic, lam**3)), odd_seed)
    return Field(even, even_correction, odd, beta, odd_cubic, odd_seed, length, height, epsilon)
