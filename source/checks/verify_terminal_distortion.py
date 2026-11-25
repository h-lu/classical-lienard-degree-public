"""Exact tests of the contracting, logistic and oblique-exit passage.

A=1, q=1/4, epsilon=1/(k^2 log 2), delta=1/k. The common oblique-exit
factor cancels. All exponentials in the remaining formula are integer
powers of two. Compare rationals without evaluating a logarithm.
"""

from fractions import Fraction


def power_of_two(exponent):
    return Fraction(2)**exponent


def repelling_exit(distance, growth):
    return growth * distance / (1 + (growth - 1) * distance)


def integrate_reciprocal(distance, steps):
    # v=1/x solves v' = -(D/epsilon)*(v-1). Each time slice halves
    # v-1. This evaluates the orbit without using the logistic quotient.
    v = 1 / distance
    for _ in range(steps):
        v = (v + 1) / 2
    return 1 / v


def main():
    q = Fraction(1, 4)
    for k in (8, 16, 32, 64):
        n = k*k
        d0 = q * power_of_two(-n)
        for budget in (Fraction(1, 2), Fraction(1), Fraction(3, 2)):
            steps = int(budget*n)
            growth = power_of_two(steps)
            baseline = repelling_exit(d0, growth)
            for y in (-1, 0, 1):
                distance = d0 * power_of_two(k*y)
                exit_value = integrate_reciprocal(distance, steps)
                assert exit_value == repelling_exit(distance, growth)
                ratio = exit_value / baseline
                limit = {Fraction(1, 2): y, Fraction(1): min(y, 0), Fraction(3, 2): 0}[budget]
                # This is |(epsilon/delta) log(ratio) - limit| < 3/k.
                assert power_of_two(k*limit - 3) < ratio < power_of_two(k*limit + 3)
            # Exact derivative ratio K'(0)/K'(d0), including the critical case.
            # For an autonomous scalar flow away from equilibria its
            # derivative equals the ratio of the exit and entry speeds.
            derivative_at_d0 = baseline*(1-baseline) / (d0*(1-d0))
            derivative_ratio = growth / derivative_at_d0
            assert derivative_ratio == (1 + (growth - 1)*d0)**2
            if budget == 1:
                assert derivative_ratio == (1 + q*(1-power_of_two(-n)))**2
                assert Fraction(3, 2) < derivative_ratio < Fraction(25, 16)
            elif budget < 1:
                assert derivative_ratio - 1 < power_of_two(-n//2)
        print(f"k={k}: nine exact gap ratios; scaled errors < 3/{k}; critical derivative ratio in (3/2,25/16).")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()
