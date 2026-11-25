"""Build the standalone reproducible certificate archive."""

import gzip
import hashlib
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SINGLE_FILES = (
    "build.sh", "replay-all.sh", "REPLAY.txt", "REPLAY-GUIDE.md",
    "tex/existence-table.tex", "tex/cubic-table.tex", "tex/wells.tex", "tex/fixed-six-table.tex",
    "verification/check_fixed_inputs.py", "verification/check_action_coordinates.py",
    "verification/fixed-degree12-replay.txt", "verification/fixed-inputs.txt",
    "verification/action-coordinates.txt",
)


CHECKS = (
    "verify_entrance.py",
    "verify_fixed_sextic.py",
    "verify_joint_shoulder_cap.py",
    "verify_one_sided_wells.py",
    "verify_outer.py",
    "verify_paper_tables.py",
    "verify_primary_sources.py",
    "verify_remote_odd_powers.py",
    "verify_seed.py",
    "verify_square_shoulders.py",
    "verify_terminal_distortion.py",
    "verify_unequal_cubic_layers.py",
    "verify_weak_return.py",
)


def scientific_inputs():
    paths = [ROOT / name for name in SINGLE_FILES]
    paths.extend(ROOT / "checks" / name for name in CHECKS)
    fixed = ROOT / "verification/fixed-degree12"
    for name in (
        "PROOF.md", "replay.sh", "compact.py", "field.py", "fixed_outer.py",
        "outer.py", "rational.py", "shooting.py", "reference/checks/verify_seed.py",
        "inputs/compact-intervals.json", "inputs/fixed-diagnostic-field.json",
        "inputs/outer-variant-field.json",
    ):
        paths.append(fixed / name)
    return sorted(paths)


CHECKSUM = "CERTIFICATE-SHA256SUMS"


def certificate_paths():
    return sorted([*scientific_inputs(), ROOT / CHECKSUM])


def main():
    manifest = "".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT).as_posix()}\n"
        for path in scientific_inputs()
    )
    (ROOT / CHECKSUM).write_text(manifest)
    with (ROOT / "certificates.tar.gz").open("wb") as output:
        with gzip.GzipFile(fileobj=output, filename="", mode="wb", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT) as tar:
                for path in certificate_paths():
                    name = path.relative_to(ROOT).as_posix()
                    info = tarfile.TarInfo(f"certificates/{name}")
                    info.size = path.stat().st_size
                    info.mode = 0o755 if name.endswith(".sh") else 0o644
                    with path.open("rb") as content:
                        tar.addfile(info, content)


if __name__ == "__main__":
    main()
