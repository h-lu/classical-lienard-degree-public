"""Check the printed comparison rows and the degree bookkeeping used to read them."""

from fractions import Fraction as Q
from pathlib import Path
import json
import re
import sys


def comparison_rows():
    degrees = range(5, 14)
    published = [2 if n == 5 else n - 2 for n in degrees]
    isolated = [None if n == 5 else (n - 1) // 2 + 4 for n in degrees]
    derived = [max(p, q or 0) for p, q in zip(published, isolated)]
    # Counts proved in S1, including the special joint construction at twelve.
    attained = [2, 6, 7, 9, 10, 14, 15, 19, 22]
    return [published, isolated, derived, attained,
            [a - b for a, b in zip(attained, derived)]]


def increasing_root(fn, lo, hi):
    lo, hi = Q(lo), Q(hi)
    assert fn(lo) <= 0 <= fn(hi)
    for _ in range(90):
        mid = (lo + hi) / 2
        if fn(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def cubic_upper(a, shift, level, witness=None):
    # a = -3 k^2. At the landing, P(-2k) = -138 k^5/5,
    # and the deleted jump has action 162 k^5/5.
    k2 = Q(-a, 3)
    kl, kh = increasing_root(lambda k: k*k - k2, Q(1, 2), 1)
    if witness is not None:
        kl, kh = ((Q(57735, 10**5), Q(57736, 10**5)) if a == -1 else (Q(1), Q(1)))
        assert kl**2 <= k2 <= kh**2
    p = lambda v: 9*v**5/5 + 2*a*v**3 + a*a*v
    offset = level - shift
    if offset**2 < 4*k2**3:
        assert witness in (None, "L")
        return -138*kh**5/5, -138*kl**5/5
    assert offset**2 > 4*k2**3
    above = offset > 0
    lo, hi = (kh, 128) if above else (-128, -kh)
    x, y = increasing_root(lambda v: v**3 + a*v - offset, lo, hi)
    if witness is not None:
        j = int(witness)
        x, y = Q(j, 10**5), Q(j+1, 10**5)
        assert (above and x > kh) or (not above and y < -kh)
        assert x**3+a*x < offset < y**3+a*y
    if above:
        return p(x) - 162*kh**5/5, p(y) - 162*kl**5/5
    return p(x), p(y)


def cubic_table(source):
    def row(label):
        line = next(line for line in source.splitlines() if line.startswith(label + " &"))
        return [Q(x.strip().strip("$").replace("10^5", "100000"))
                for x in line.removesuffix(r"\\").split("&")[1:]]
    levels, lower, upper = row(r"$\xi$"), row("Lower bound"), row("Upper bound")
    assert len(levels) == len(lower) == len(upper) == 8
    def witnesses(label):
        line = next(line for line in source.splitlines() if line.startswith(label + " &"))
        return [c.strip().strip("$") for c in line.removesuffix(r"\\").split("&")[1:]]
    plus, minus = witnesses("$j_+$"), witnesses("$j_-$")
    assert len(plus) == len(minus) == len(levels)
    for t, a, b, jp, jm in zip(levels, lower, upper, plus, minus):
        for wp, wm in ((None, None), (jp, jm)):
            pl, ph = cubic_upper(-1, -8, t, wp)
            ml, mh = cubic_upper(-3, 8, t, wm)
            assert a < pl - mh - Q(493, 2) <= ph - ml - Q(493, 2) < b
    print("Independent inversion and the printed root brackets give all eight split-cubic enclosures.")


def well_table(source):
    line = next(line for line in source.splitlines() if line.strip().startswith("100B_d &"))
    bands = [(int(a), int(b)) for a, b in re.findall(r"\((-?\d+),(-?\d+)\)", line)]
    def numerators(key):
        line = next(line for line in source.splitlines() if line.startswith(key+" &"))
        return [int(c) for c in line.removesuffix(r"\\").split("&")[1:]]
    xs, vs = numerators("X"), numerators("V")
    assert len(bands) == len(xs) == len(vs) == 7
    for d, ((a, b), nx, nv) in enumerate(zip(bands, xs, vs), 3):
        m = 2*d
        t = Q(d, d+1)**d / (d+1)
        xl, xh = increasing_root(lambda x: x**d*(1+x)-t, 0, 1)
        vl, vh = increasing_root(lambda v: v**d*(v-1)-t, 1, 2)

        def action(s, start, end):
            coeffs = (m, 2*s*(m+1), m+2)
            total = Q(0)
            for i, ci in enumerate(coeffs):
                for j, cj in enumerate(coeffs):
                    power = 2*m - 2 + i + j
                    total += ci*cj*(end**power - start**power)/power
            return total

        n = 10**6
        assert nx**d*(n+nx) < n**(d+1)*t < (nx+1)**d*(n+nx+1)
        assert nv**d*(nv-n) < n**(d+1)*t < (nv+1)**d*(nv+1-n)
        for x0, x1, v0, v1 in ((xl, xh, vl, vh),
                               (Q(nx,n), Q(nx+1,n), Q(nv,n), Q(nv+1,n))):
            low = 100*(1-action(-1, 1, v1)/action(1, 0, x0))
            high = 100*(1-action(-1, 1, v0)/action(1, 0, x1))
            assert a < low <= high < b
    print("Printed rational root brackets and independent integration give the seven well bands (m=6 through 18).")


def fixed_table(root):
    rows = (root/"tex/fixed-six-table.tex").read_text().split(r"\midrule")[1].split(r"\bottomrule")[0]
    fixed = root/"verification"
    intervals = json.loads((fixed/"fixed-degree12/inputs/compact-intervals.json").read_text())
    output = (fixed/"fixed-degree12-replay.txt").read_text().split("Fixed rational variant with theta=-1e-20")[1]
    results = re.findall(r"^  \(([\d.]+), ([\d.]+)\): signs ([+-])1/([+-])1, integral difference in \(([^,]+), ([^)]+)\)", output, re.MULTILINE)
    lines = [line for line in rows.splitlines() if "&" in line]
    assert len(results) == len(intervals) == len(lines) == 5
    for line, (l, r), (lo, hi, left, right, a, b) in zip(lines, intervals, results):
        h0, h1, s0, s1, bounds = [cell.strip().strip("$") for cell in line.removesuffix(r"\\").split("&")]
        assert (Q(h0), Q(h1)) == (Q(lo), Q(hi)) == (Q(l), Q(r))
        assert (s0, s1) == (left, right)
        assert tuple(map(Q,bounds.strip("()").split(","))) == (10**5*Q(a), 10**5*Q(b))
    print("The five printed compact brackets, endpoint signs and full-interval divergence bounds match the fixed-field certificate.")


def main():
    table = Path(__file__).resolve().parents[1] / "tex" / "existence-table.tex"
    text = table.read_text().split(r"\midrule")[1].split(r"\bottomrule")[0]
    lines = [line for line in text.splitlines() if "&" in line]
    cubic_table((table.parent / "cubic-table.tex").read_text())
    well_table((table.parent / "wells.tex").read_text())
    fixed_table(table.parent.parent)
    expected = comparison_rows()
    proof_cells = lines.pop().removesuffix(r"\\").split("&")[1:]
    evidence_cells = lines.pop().removesuffix(r"\\").split("&")[1:]
    statement_cells = lines.pop().removesuffix(r"\\").split("&")[1:]
    octic = r"\ref{cor:octic-low}"
    assert statement_cells == [r"\hyperref[thm:weak]{A.1(iii)}", "S1.2", "S1.2", octic,
                               "S1.2", "S1.2", octic, "S1.1", octic]
    assert evidence_cells == ["A", "V", "V", "A", "V", "V", "A", "V", "A"]
    assert proof_cells == [
        "$d=2$", "S3.1", "S3.1", r"\hyperlink{octic-core}{core}",
        "S5.2", "S4", r"\hyperlink{odd-closure}{odd, $t=0$}", "S5.1",
        r"\hyperlink{odd-closure}{odd, $t=1$}",
    ]
    assert len(lines) == len(expected)
    for line, row in zip(lines, expected):
        label, *cells = line.removesuffix(r"\\").split("&")
        values = [None if cell.strip() == "--" else int(cell) for cell in cells]
        assert values == row, label
        print(f"{label.strip()}: " + " ".join("--" if n is None else str(n) for n in values))
    # The even closure starts from eighteen. Its first well therefore gives 25.
    for n, count in ((13, 22), (14, 25), (15, 29), (16, 32), (17, 36)):
        odd = n % 2
        t = (n - (11 if odd else 12)) // 2
        assert (15 if odd else 18) + 7*t == count == (7*n - 47)//2
    print("Octic closures: 11/15, 12/18; one well: 13/22, 14/25; next: 15/29, 16/32, 17/36.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run Python with assertions enabled.")
    main()
