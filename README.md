# Classical polynomial Liénard systems at fixed degree

This repository is the reproducibility companion to Haibo Lu's preprint, [*Relative return gaps and a 7/2 asymptotic lower bound for hyperbolic limit cycles in classical Liénard systems*](https://doi.org/10.5281/zenodo.24517779). The two PDFs and the downloadable source archive match the files deposited on Zenodo.

The [article](main.pdf) proves a lower bound of \(\lceil(7n-48)/2\rceil\) coexisting hyperbolic limit cycles at each exact primitive degree \(n\ge14\). The [electronic supplement](low-degree-supplement.pdf) gives the low-degree results, including a nineteen-cycle existence theorem at degree twelve and a six-cycle certificate for a specified rational field. A second rational field attains five cycles in the weak full-return domain.

## Reproduce the article and supplement

The complete editable archive is `source-reproduction.tar.gz`. Place the accompanying `.sha256` file beside it and run the following commands.

```sh
sha256sum -c source-reproduction.tar.gz.sha256
tar -xzf source-reproduction.tar.gz
cd classical-lienard-degree
./reproduce.sh
```

The driver checks the scientific inventory, rebuilds both PDFs from the [LaTeX sources](source/main.tex), regenerates the [certificate archive](source/certificates.tar.gz) and replays its four reference transcripts. The [build instructions](source/BUILD.md) list the dependencies and single-document commands. For the finite checks alone, extract the certificate archive and follow its [replay guide](source/REPLAY-GUIDE.md).

The [certificate index](CERTIFICATE-INDEX.md) links the exact inputs and outputs to their mathematical use. The [mechanism note](MECHANISM.md) gives the relative-gap estimate, its hypotheses, degree cost and a scalar passage illustrating its limit of applicability.

## Citation and reuse

Lu, Haibo (2026). *Relative return gaps and a 7/2 asymptotic lower bound for hyperbolic limit cycles in classical Liénard systems*. Zenodo. [doi:10.5281/zenodo.24517779](https://doi.org/10.5281/zenodo.24517779).

Use the Zenodo DOI to cite the preprint, and include GitHub tag `v1.0.0` if referring to this exact mirror. The article, supplement, original certificate inputs and outputs, and original reproducibility source mirror the Zenodo package under [CC BY 4.0](LICENSE). External dependencies retain their own licenses.
