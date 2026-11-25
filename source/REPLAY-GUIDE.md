# Standalone certificate replay

Extract the separately supplied archive in a fresh directory.

```sh
tar -xzf certificates.tar.gz
cd certificates
sha256sum -c CERTIFICATE-SHA256SUMS
./replay-all.sh
python3 -B verification/check_action_coordinates.py > coordinates.actual
cmp verification/action-coordinates.txt coordinates.actual
```

Bash, Python 3.9 or later with its standard library, tar, cmp and sha256sum suffice. The checksum verifies the expanded payload against the enclosed inventory. The replay compares all thirteen manuscript programs, four fixed-field programs and the coefficient and escape test against three expected transcripts. The final two commands check the action-coordinate identities using both rational fields and compare the fourth transcript. Supplement Appendix SB proves the six-cycle result and derives the fixed-scale matching estimates. The certificate archive also contains their derivations in [PROOF.md](verification/fixed-degree12/PROOF.md). The nineteen-cycle result uses validated sextic data and analytic matching. The four fixed-field programs establish six hyperbolic cycles for the specified negative-tilt polynomial.
