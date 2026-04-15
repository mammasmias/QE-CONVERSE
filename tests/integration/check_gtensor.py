#!/usr/bin/env python3
"""
Check EPR g-tensor outputs against reference values.

Usage:
    check_gtensor.py [options] output1.out [output2.out ...]
    check_gtensor.py [options] --outdir DIR

For each output file the script looks for a file with the same name in the
reference directory and compares all g-tensor quantities within tolerances.

Options:
    --outdir DIR     Directory containing output *.out files to check
    --refdir DIR     Reference directory [default: reference/ next to outputs]
    --atol-rmc X     Absolute tolerance for RMC quantities (ppm) [default: 0.05]
    --atol-so  X     Absolute tolerance for SO/total quantities (ppm) [default: 1.0]
    --atol-mag X     Absolute tolerance for magnetization (a.u.) [default: 1e-4]
    --rtol X         Relative tolerance applied to all quantities [default: 1e-4]

Exit code: 0 if all checks pass, 1 if any fail.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from parse_output import parse_gtensor, parse_metadata, parse_orbital_magnetization

# Default tolerances
ATOL_RMC    = 0.05   # ppm
ATOL_SO     = 1.0    # ppm
ATOL_MAG    = 1e-4   # a.u.
ATOL_ENERGY = 1e-5   # Ry
RTOL        = 1e-4


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
                f"       got {val:.6f}  ref {ref:.6f}"
                f"  |diff|={abs(val - ref):.4g}  atol={atol}")
    return True, f"  OK   {name}"


def check_output(out_path, ref_path, *,
                 atol_rmc=ATOL_RMC, atol_so=ATOL_SO,
                 atol_mag=ATOL_MAG, atol_energy=ATOL_ENERGY, rtol=RTOL):
    """
    Compare one output file against its reference.

    Returns (passed: bool, messages: list[str]).
    """
    out_text = out_path.read_text()
    ref_text = ref_path.read_text()

    result   = parse_gtensor(out_text)
    ref      = parse_gtensor(ref_text)
    meta     = parse_metadata(out_text)
    ref_meta = parse_metadata(ref_text)
    orb      = parse_orbital_magnetization(out_text)
    orb_ref  = parse_orbital_magnetization(ref_text)

    messages = []
    all_ok = True

    # Convergence
    if not meta.get("converged"):
        messages.append("  FAIL convergence: JOB DONE not found in output")
        all_ok = False
    else:
        messages.append("  OK   convergence")

    # SCF total energy
    ok, msg = _check_close("total energy (Ry)",
                           meta.get("total_energy"), ref_meta.get("total_energy"),
                           atol=atol_energy, rtol=rtol)
    messages.append(msg)
    if not ok:
        all_ok = False

    # Intermediate magnetization contributions (without Berry curvature)
    for key in ("delta_m_bare", "delta_m_para", "delta_m_dia",
                "m_lc_no_bc", "m_ic_no_bc", "m_tot_no_bc"):
        ok, msg = _check_close(key, orb.get(key), orb_ref.get(key),
                               atol=atol_mag, rtol=rtol)
        messages.append(msg)
        if not ok:
            all_ok = False

    # Final g-tensor quantities (use parse_gtensor values which include
    # m_lc/m_ic/m_tot with Berry curvature via last-occurrence matching)
    checks = [
        ("delta_g_rmc",       atol_rmc),
        ("delta_g_rmc_gipaw", atol_rmc),
        ("delta_g_so",        atol_so),
        ("delta_g_total",     atol_so),
        ("m_lc",              atol_mag),
        ("m_ic",              atol_mag),
        ("m_tot",             atol_mag),
    ]
    for key, atol in checks:
        ok, msg = _check_close(key, result.get(key), ref.get(key),
                               atol=atol, rtol=rtol)
        messages.append(msg)
        if not ok:
            all_ok = False

    return all_ok, messages


def main():
    parser = argparse.ArgumentParser(
        description="Check EPR g-tensor outputs against reference values.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--files", nargs="+", required=True,
                        help="Output file(s) to check")
    parser.add_argument("--refdir",
                        help="Reference directory [default: reference/ next to outputs]")
    parser.add_argument("--atol-rmc", type=float, default=ATOL_RMC,
                        dest="atol_rmc",
                        help=f"Tolerance for RMC quantities (ppm) [default: {ATOL_RMC}]")
    parser.add_argument("--atol-so", type=float, default=ATOL_SO,
                        dest="atol_so",
                        help=f"Tolerance for SO/total quantities (ppm) [default: {ATOL_SO}]")
    parser.add_argument("--atol-mag", type=float, default=ATOL_MAG,
                        dest="atol_mag",
                        help=f"Tolerance for magnetization (a.u.) [default: {ATOL_MAG}]")
    parser.add_argument("--atol-energy", type=float, default=ATOL_ENERGY,
                        dest="atol_energy",
                        help=f"Tolerance for SCF total energy (Ry) [default: {ATOL_ENERGY}]")
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

    for out_path in out_files:
        ref_path = ref_dir / out_path.name
        if not ref_path.exists():
            print(f"\n[SKIP] {out_path.name}: no reference file at {ref_path}")
            n_skip += 1
            continue

        print(f"\n[CHECK] {out_path.name}")
        passed, messages = check_output(
            out_path, ref_path,
            atol_rmc=args.atol_rmc, atol_so=args.atol_so,
            atol_mag=args.atol_mag, atol_energy=args.atol_energy,
            rtol=args.rtol)
        for msg in messages:
            print(msg)
        if passed:
            n_pass += 1
            print("  → PASS")
        else:
            n_fail += 1
            print("  → FAIL")

    print(f"\nSummary: {n_pass} passed, {n_fail} failed, {n_skip} skipped")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
