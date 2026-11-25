# Building the article and its electronic supplement

From the extracted `classical-lienard-degree/` directory, the following command verifies [the inventory](../SHA256SUMS), builds both documents, regenerates the numerical archive and replays its certificates. It compares both PDFs and the archive with the distributed copies and compares the new outputs with the four reference transcripts.

```sh
./reproduce.sh
```

The build uses Bash, Python 3.9 or later (standard library), tar, cmp, sha256sum, pdfLaTeX and BibTeX. The TeX inputs require Latin Modern, the AMS packages, mathtools, microtype, geometry, TikZ, booktabs, array and hyperref. TeX Live 2025/Debian gives the distributed PDF bytes with these sources. The build fixes the epoch at 1795392000, uses UTC, disables shell escape and checks the final typesetting and bibliography logs. Other TeX releases can give different PDF bytes.

For an editable working copy, run `./source/build.sh main` or `./source/build.sh supplement`. This writes the selected PDF alongside its root TeX file in `source/`. Both documents have their own bibliography and internal references.

To run the scientific checks only, see the [certificate replay guide](REPLAY-GUIDE.md). The [certificate index](../CERTIFICATE-INDEX.md) identifies each calculation and its role in the proof.
