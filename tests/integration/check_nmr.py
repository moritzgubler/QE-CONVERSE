#!/usr/bin/env python3
"""
Check NMR chemical shift outputs against reference values.

Each QE-CONVERSE NMR run computes the chemical shift for one atom along one
Cartesian direction (m_0(dir) = 1.0).  The output contains a single
"Chemical shift (ppm):" line with three components.

This script compares each output file against the corresponding reference file
(matched by filename) and optionally computes the isotropic chemical shift
from groups of x/y/z direction runs.

Usage:
    check_nmr.py [options] output1.out [output2.out ...]
    check_nmr.py [options] --outdir DIR

Options:
    --outdir DIR      Directory containing output *.out files to check
    --refdir DIR      Reference directory [default: reference/ next to outputs]
    --atol-shift X    Absolute tolerance for shift components (ppm) [default: 2.0]
    --atol-iso X      Absolute tolerance for isotropic shift (ppm) [default: 1.0]
    --rtol X          Relative tolerance applied to all quantities [default: 1e-3]

Exit code: 0 if all checks pass, 1 if any fail.
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from parse_output import parse_nmr_shift, parse_metadata

# Default tolerances
ATOL_SHIFT = 2.0   # ppm, individual components
ATOL_ISO   = 1.0   # ppm, isotropic
RTOL       = 1e-3


def _check_close(name, val, ref, atol, rtol=RTOL):
    """Return (ok, message).  val/ref may be scalars or equal-length lists."""
    if val is None:
        return False, f"  FAIL {name}: could not parse value from output"
    if ref is None:
        return False, f"  FAIL {name}: could not parse value from reference"

    if isinstance(val, list):
        pairs = list(zip(val, ref))
        ok = all(abs(v - r) <= atol + rtol * abs(r) for v, r in pairs)
        if not ok:
            diffs = [f"{abs(v - r):.4g}" for v, r in pairs]
            return False, (
                f"  FAIL {name}:\n"
                f"       got {[f'{v:.4f}' for v in val]}\n"
                f"       ref {[f'{r:.4f}' for r in ref]}\n"
                f"       |diff| {diffs}  atol={atol}")
    else:
        ok = abs(val - ref) <= atol + rtol * abs(ref)
        if not ok:
            return False, (
                f"  FAIL {name}:\n"
                f"       got {val:.4f}  ref {ref:.4f}"
                f"  |diff|={abs(val - ref):.4g}  atol={atol}")
    return True, f"  OK   {name}"


def check_output(out_path, ref_path, *, atol_shift=ATOL_SHIFT, rtol=RTOL):
    """
    Compare one output file against its reference.

    Returns (passed: bool, messages: list[str]).
    """
    out_text = out_path.read_text()
    ref_text = ref_path.read_text()

    result = parse_nmr_shift(out_text)
    ref    = parse_nmr_shift(ref_text)
    meta   = parse_metadata(out_text)

    messages = []
    all_ok = True

    # Convergence
    if not meta.get("converged"):
        messages.append("  FAIL convergence: JOB DONE not found in output")
        all_ok = False
    else:
        messages.append("  OK   convergence")

    # Chemical shift vector
    ok, msg = _check_close("shift (ppm)", result.get("shift"), ref.get("shift"),
                           atol=atol_shift, rtol=rtol)
    messages.append(msg)
    if not ok:
        all_ok = False

    return all_ok, messages


def _group_by_atom(file_paths):
    """
    Group output paths by atom label (e.g. 'Si1', 'O4').

    Expects filenames matching *_<label><dir>.out where <dir> is x, y, or z.
    Returns dict: label → {'x': path, 'y': path, 'z': path}.
    """
    groups = {}
    for p in file_paths:
        m = re.match(r'.*?([A-Za-z]+\d+)([xyz])\.out$', p.name)
        if m:
            label, direction = m.group(1), m.group(2)
            groups.setdefault(label, {})[direction] = p
    return groups


def compute_isotropic(x_text, y_text, z_text):
    """
    Compute isotropic shift from three single-direction outputs.

    The isotropic shielding is (sigma_xx + sigma_yy + sigma_zz) / 3,
    where sigma_ii is the i-th component from the i-direction run.
    The core shift (constant per element) is added if present.

    Returns (sigma_iso_bare, sigma_iso_total) in ppm, or (None, None).
    """
    rx = parse_nmr_shift(x_text)
    ry = parse_nmr_shift(y_text)
    rz = parse_nmr_shift(z_text)

    sx = rx.get("shift")
    sy = ry.get("shift")
    sz = rz.get("shift")
    if sx is None or sy is None or sz is None:
        return None, None

    # Diagonal elements: first component from x-run, second from y-run, third from z-run
    sigma_iso = (sx[0] + sy[1] + sz[2]) / 3.0

    core = rx.get("core") or ry.get("core") or rz.get("core")
    sigma_total = sigma_iso + core if core is not None else None

    return sigma_iso, sigma_total


def main():
    parser = argparse.ArgumentParser(
        description="Check NMR chemical shift outputs against reference values.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--files", nargs="+", required=True,
                        help="Output file(s) to check")
    parser.add_argument("--refdir",
                        help="Reference directory [default: reference/ next to outputs]")
    parser.add_argument("--atol-shift", type=float, default=ATOL_SHIFT,
                        dest="atol_shift",
                        help=f"Tolerance for shift components (ppm) [default: {ATOL_SHIFT}]")
    parser.add_argument("--atol-iso", type=float, default=ATOL_ISO,
                        dest="atol_iso",
                        help=f"Tolerance for isotropic shift (ppm) [default: {ATOL_ISO}]")
    parser.add_argument("--rtol", type=float, default=RTOL,
                        help=f"Relative tolerance [default: {RTOL}]")
    args = parser.parse_args()

    out_files = [Path(f) for f in args.files]

    # Determine reference directory
    if args.refdir:
        ref_dir = Path(args.refdir)
    else:
        ref_dir = out_files[0].parent / "reference"
    if not ref_dir.is_dir():
        print(f"ERROR: reference directory not found: {ref_dir}",
              file=sys.stderr)
        sys.exit(1)

    n_pass = n_fail = n_skip = 0

    # Per-file component checks
    for out_path in out_files:
        ref_path = ref_dir / out_path.name
        if not ref_path.exists():
            print(f"\n[SKIP] {out_path.name}: no reference file at {ref_path}")
            n_skip += 1
            continue

        print(f"\n[CHECK] {out_path.name}")
        passed, messages = check_output(
            out_path, ref_path,
            atol_shift=args.atol_shift, rtol=args.rtol)
        for msg in messages:
            print(msg)
        if passed:
            n_pass += 1
            print("  → PASS")
        else:
            n_fail += 1
            print("  → FAIL")

    # Isotropic shift summary (groups x/y/z runs per atom)
    out_groups = _group_by_atom(out_files)
    ref_groups = _group_by_atom(sorted(ref_dir.glob("*.out")))

    if out_groups:
        print("\n--- Isotropic chemical shifts ---")
        for label in sorted(out_groups):
            og = out_groups[label]
            rg = ref_groups.get(label, {})
            if not all(d in og for d in ("x", "y", "z")):
                print(f"  {label}: skipping isotropic (not all directions present)")
                continue

            iso_out, tot_out = compute_isotropic(
                og["x"].read_text(), og["y"].read_text(), og["z"].read_text())

            ref_iso = ref_tot = None
            if all(d in rg for d in ("x", "y", "z")):
                ref_iso, ref_tot = compute_isotropic(
                    rg["x"].read_text(), rg["y"].read_text(), rg["z"].read_text())

            if iso_out is not None and tot_out is not None:
                print(f"  {label}: sigma_bare={iso_out:.4f} ppm, "
                      f"sigma_total={tot_out:.4f} ppm", end="")
                if ref_tot is not None:
                    ok, msg = _check_close(
                        f"{label} isotropic", tot_out, ref_tot,
                        atol=args.atol_iso, rtol=args.rtol)
                    if ok:
                        print(f"  (ref={ref_tot:.4f}, OK)")
                    else:
                        print(f"  (ref={ref_tot:.4f}, FAIL diff={abs(tot_out - ref_tot):.4g})")
                        n_fail += 1
                        n_pass -= 1  # adjust since file-level already counted
                else:
                    print()
            else:
                print(f"  {label}: could not compute isotropic shift")

    print(f"\nSummary: {n_pass} passed, {n_fail} failed, {n_skip} skipped")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
