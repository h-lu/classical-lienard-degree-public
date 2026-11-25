"""Exact rational polynomials and outward-rounded dyadic intervals."""

from fractions import Fraction as Q
from math import comb


class Box:
    scale = 2**4096

    def __init__(self, lower=0, upper=None):
        if upper is None:
            upper = lower
        self._lower = (Q(lower) * self.scale).__floor__()
        self._upper = (Q(upper) * self.scale).__ceil__()
        assert self._lower <= self._upper

    @classmethod
    def _from_grid(cls, lower, upper):
        new = cls.__new__(cls)
        new._lower = lower
        new._upper = upper
        assert lower <= upper
        return new

    @property
    def lower(self):
        return Q(self._lower, self.scale)

    @property
    def upper(self):
        return Q(self._upper, self.scale)

    @staticmethod
    def wrap(value):
        return value if isinstance(value, Box) else Box(value)

    def __neg__(self):
        return self._from_grid(-self._upper, -self._lower)

    def __add__(self, other):
        other = self.wrap(other)
        return self._from_grid(self._lower + other._lower, self._upper + other._upper)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-self.wrap(other))

    def __rsub__(self, other):
        return self.wrap(other) + (-self)

    def __mul__(self, other):
        other = self.wrap(other)
        corners = (self._lower * other._lower, self._lower * other._upper,
                   self._upper * other._lower, self._upper * other._upper)
        return self._from_grid(min(corners) // self.scale, -(-max(corners) // self.scale))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.wrap(other)
        assert other._lower * other._upper > 0, "division by an interval containing zero"
        square = self.scale**2
        reciprocal = self._from_grid(square // other._upper, -(-square // other._lower))
        return self * reciprocal

    def __rtruediv__(self, other):
        return self.wrap(other) / self

    def __pow__(self, power):
        assert isinstance(power, int) and power >= 0
        if power == 0:
            return Box(1)
        half = self ** (power // 2)
        even = half * half
        return even * self if power % 2 else even

    def radius_bound(self):
        return Q(max(abs(self._lower), abs(self._upper)), self.scale)

    def strictly_inside(self, lower, upper):
        return Q(lower) < self.lower <= self.upper < Q(upper)


def add(left, right):
    result = [0] * max(len(left), len(right))
    for term in (left, right):
        for index, coefficient in enumerate(term):
            result[index] += coefficient
    return result


def scale(poly, factor):
    return [factor * coefficient for coefficient in poly]


def sub(left, right):
    return add(left, scale(right, -1))


def mul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def deriv(poly):
    return [i * coefficient for i, coefficient in enumerate(poly) if i]


def primitive(poly):
    return [0, *(Q(1, i + 1) * coefficient for i, coefficient in enumerate(poly))]


def at(poly, x):
    result = 0
    for coefficient in reversed(poly):
        result = result * x + coefficient
    return result


def norm(poly, radius):
    return sum(Box.wrap(coefficient).radius_bound() * radius**i
               for i, coefficient in enumerate(poly))


def translate(poly, x):
    return [sum(poly[j] * comb(j, i) * x**(j - i) for j in range(i, len(poly)))
            for i in range(len(poly))]


def bernstein_positive(poly, left, right):
    degree = len(poly) - 1
    scaled = [c * (right - left)**i for i, c in enumerate(translate(poly, left))]
    return all(Box.wrap(sum(scaled[i] * Q(comb(k, i), comb(degree, i))
                            for i in range(k + 1))).lower > 0
               for k in range(degree + 1))


def positive_on(poly, left, right, cells=16):
    grid = [left + (right - left) * Q(i, cells) for i in range(cells + 1)]
    assert all(bernstein_positive(poly, a, b) for a, b in zip(grid, grid[1:])), (left, right)


def decreasing_zero(residual, left, right, steps=560):
    """Isolate every zero in a family of strictly decreasing functions."""
    assert residual(left).lower > 0 and residual(right).upper < 0
    a, b, c, d = left, right, left, right
    for _ in range(steps):
        x = (a + b) / 2
        if residual(x).lower > 0:
            a = x
        else:
            b = x
        y = (c + d) / 2
        if residual(y).upper < 0:
            d = y
        else:
            c = y
    return Box(a, d)
