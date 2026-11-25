"""Independent half-return integration at the rational centre of Palomar, App. A."""

from fractions import Fraction as Q

from verify_primary_sources import ABSOLUTE_GAPS, PALOMAR_COEFFICIENTS
from verify_seed import return_height


HEIGHTS = tuple(map(Q, (
    "1/2",
    "267947284076040944055904337196005393766617/72697234874773572860415681250000000000000",
    "70494125095724733767180560718998241791867/18174308718693393215103920312500000000000",
    "4684339993646143644823568052538315255527353/1163155757996377165766650900000000000000000",
    "2535575176391823260968901707271760975699691/581577878998188582883325450000000000000000",
    "1324047472268609322871629305970037251288413/290788939499094291441662725000000000000000",
    "2723570269909365404647409649413501806826697/581577878998188582883325450000000000000000",
)))


def main():
    coefficients = (Q(0), *PALOMAR_COEFFICIENTS)
    reflected = tuple((-1)**j * c for j, c in enumerate(coefficients))
    assert all(a < b for a, b in zip(HEIGHTS, HEIGHTS[1:]))
    for i, (height, margin) in enumerate(zip(HEIGHTS, ABSOLUTE_GAPS)):
        plus, psteps = return_height(coefficients, Q(0), height)
        minus, msteps = return_height(reflected, Q(0), height)
        lower = plus[0] - minus[1] if i % 2 == 0 else minus[0] - plus[1]
        assert lower > Q(margin), (i, lower)
        sign = "+" if i % 2 == 0 else "-"
        print(f"h{i}: {sign}, signed magnitude > {margin}, {psteps} + {msteps} Taylor steps",
              flush=True)
    print("Fixed rational sextic: six disjoint brackets for isolated cycles.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()
