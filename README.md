# Relative return gaps in classical Liénard systems

**Haibo Lu**, Shanghai Institute of Technology. [Zenodo preprint and permanent archive](https://doi.org/10.5281/zenodo.24517779) · [GitHub release v1.0.0](https://github.com/h-lu/classical-lienard-degree-public/releases/tag/v1.0.0).

This repository accompanies the preprint *Relative return gaps and a 7/2 asymptotic lower bound for hyperbolic limit cycles in classical Liénard systems*. It contains the signed article, its electronic supplement, editable LaTeX, and the exact inputs and programs needed to replay the finite checks.

## Problem and main result

The paper studies the classical Liénard equation with linear restoring force, written as the planar system

```math
\dot{x}=y-F(x),\qquad \dot{y}=-x,\qquad F(x)=\sum_{j=1}^{n}a_jx^j,\quad a_n\ne0.
```

Here **n is the exact degree of the primitive F**. The damping polynomial F′ has degree n−1. The budget is polynomial degree, not the number of nonzero terms. A hyperbolic cycle has a transverse return multiplier different from one.

**Theorem 1.1.** For each integer n ≥ 14 there is a nonempty open set of exact-degree coefficient vectors such that every field in that set has at least

```math
L(n)=\left\lceil\frac{7n-48}{2}\right\rceil
```

geometrically distinct hyperbolic limit cycles in the **same** field. The open set contains rational coefficient vectors. This is a constructive lower bound, with asymptotic coefficient 7/2; it does not assert a maximum number of cycles.

### Low degrees and the degree-twelve results

The following are the lower bounds proved in the article and its supplement. See Table 1 of the article for the separate published and preprint comparisons, proof locations, and the type of evidence used in each column.

| Exact degree of F | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Coexisting hyperbolic cycles, at least | 2 | 6 | 7 | 9 | 10 | 14 | 15 | 19 | 22 |

- **Corollary 1.2** isolates the octic construction: it yields 18 cycles at degree 12. **Supplement Theorem S1.1** uses a separate construction and raises the degree-12 lower bound to 19.
- **The 19-cycle theorem** combines finite enclosures with analytic matching, simultaneous parameter choice, and hyperbolicity. It is an existence result on an open degree-12 coefficient set.
- **The fixed rational certificates are separate.** Supplement Proposition SB.1 proves at least six hyperbolic cycles for a specified rational field. Supplement Remark SA.1 gives five in the weak full-return domain for a different specified field.

## What the construction uses

The proof compares return heights after reflecting one half-plane and reversing time there. Section 3 gives a relative-gap estimate; Sections 4–5 use shallow wells and higher-degree terms to create new sign changes and close the corresponding returns on separated scales. The finite checks verify such inputs as sign tables, shooting and entrance enclosures, and fixed-field divergence bounds. They complement the analytic arguments that match the parameters and turn the signs into simultaneous hyperbolic cycles.

## Read, download, and reproduce

| File or guide | Start here for |
| --- | --- |
| [Article PDF](main.pdf) | Theorem 1.1, Corollary 1.2, literature comparison, and the construction |
| [Electronic supplement PDF](low-degree-supplement.pdf) | Low-degree proofs and the specified rational fields |
| [Mechanism note](MECHANISM.md) | Hypotheses, degree cost, and limits of the relative-gap estimate |
| [Certificate index](CERTIFICATE-INDEX.md) | Which exact input, program, and transcript supports each claim |
| [Editable source archive](source-reproduction.tar.gz) · [SHA-256](source-reproduction.tar.gz.sha256) | A self-contained copy of both TeX sources and the reproduction files |
| [Standalone certificate archive](source/certificates.tar.gz) · [Replay guide](source/REPLAY-GUIDE.md) | The finite checks without rebuilding the PDFs |

To check the downloaded archive, build both PDFs, and replay the certificates:

```bash
sha256sum -c source-reproduction.tar.gz.sha256
tar -xzf source-reproduction.tar.gz
cd classical-lienard-degree
./reproduce.sh
```

The same `./reproduce.sh` works at the root of a Git clone. It verifies the scientific inventory, builds the article and supplement, regenerates the certificate archive, and compares four replay transcripts. See [build instructions](source/BUILD.md) for prerequisites and the reference TeX version. Other TeX distributions can change PDF bytes without changing the mathematical output.

## Citation and reuse

Lu, Haibo (2026). *Relative return gaps and a 7/2 asymptotic lower bound for hyperbolic limit cycles in classical Liénard systems*. Zenodo. [doi:10.5281/zenodo.24517779](https://doi.org/10.5281/zenodo.24517779). Machine-readable metadata: [CITATION.cff](CITATION.cff). Cite the Zenodo DOI for the preprint and the release tag when referring to a particular GitHub snapshot.

The article, electronic supplement, original certificate inputs and outputs, and original reproducibility source mirror the public Zenodo package under [CC BY 4.0](LICENSE). External dependencies retain their own licenses.
