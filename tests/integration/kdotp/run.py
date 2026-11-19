#!/usr/bin/env python3
"""End-to-end Sternheimer regression: EPR, NMR, Hubbard and MPI pools.

Use PW and QECONVERSE for the executable paths. All calculations run in a
fresh temporary directory; source inputs and stored references are unchanged.
"""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from check_gtensor import check_output
from parse_output import parse_gtensor

PW = os.environ.get('PW', 'pw.x')
CONVERSE = os.environ.get('QECONVERSE', str(ROOT.parents[1] / 'bin/qe-converse.x'))
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')


def run(exe, inp, out, cwd, ranks=1, pools=None):
    cmd = ['mpirun', '--oversubscribe', '-np', str(ranks), exe]
    if pools:
        cmd += ['-nk', str(pools)]
    with inp.open('rb') as stdin, out.open('wb') as stdout:
        subprocess.run(cmd, cwd=cwd, stdin=stdin, stdout=stdout,
                       stderr=subprocess.STDOUT, check=True, env=ENV, timeout=1200)
    text = out.read_text()
    if 'JOB DONE' not in text or '%%%%%%%%' in text:
        raise RuntimeError(f'Incomplete calculation: {out}')


def prepare(source, dest):
    dest.mkdir()
    for inp in source.glob('*.in'):
        text = inp.read_text()
        text = re.sub(r"pseudo_dir\s*=\s*['\"]([^'\"]+)['\"]",
                      lambda m: "pseudo_dir = '" + str((source / m.group(1)).resolve()) + "'", text, flags=re.I)
        (dest / inp.name).write_text(text)


def converse_input(inp, method, step='0.005'):
    text = inp.read_text().replace("'covariant'", f"'{method}'")
    text = re.sub(r'q_gipaw\s*=\s*[\d.dDeE+-]+', f'q_gipaw = {step}', text)
    target = inp.with_name(method + '-' + inp.name)
    target.write_text(text)
    return target


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', nargs='+', choices=['epr-serial', 'epr-pools',
                        'nmr-pools', 'hubbard-smoke'], help='Run only selected case names')
    args = parser.parse_args()
    # Preserve logs on failure to make numerical regressions diagnosable.
    work = Path(tempfile.mkdtemp(prefix='qe-converse-kdotp-'))
    print(f'Test outputs: {work}', flush=True)
    serial = None
    for name, source, scf, inp, ranks, pools in [
        ('epr-serial', 'CO+', 'pw_scf.in', 'gtensor_1.in', 1, 1),
        ('epr-pools', 'CO+', 'pw_scf.in', 'gtensor_1.in', 4, 2),
        ('nmr-pools', 'nacl_nmr_kpool', 'pw.scf.in', 'Na1y.in', 8, 4),
        ('hubbard-smoke', 'licoo2_U', 'pw.scf.in', 'li1x.in', 4, 1),
    ]:
        if args.cases and name not in args.cases:
            continue
        dest = work / name
        prepare(ROOT / source, dest)
        if name == 'hubbard-smoke':
            # Cheap response regression, NOT a converged material prediction.
            scf_input = dest / scf
            text = scf_input.read_text().replace('ecutwfc = 80.0', 'ecutwfc = 30.0')
            text = text.replace('3 3 1 0 0 0', '1 1 1 1 1 1')
            scf_input.write_text(text)
        run(PW, dest / scf, dest / 'scf.out', dest, ranks)
        for method in ('covariant', 'kdotp'):
            run(CONVERSE, converse_input(dest / inp, method),
                dest / (method + '.out'), dest, ranks, pools)
        kd, cv = dest / 'kdotp.out', dest / 'covariant.out'
        residuals = re.findall(r'kdotp direction.*residual=\s*([\d.Ee+-]+)', kd.read_text())
        if not residuals or any(float(r) > 1.01e-8 for r in residuals):
            raise AssertionError(f'Missing or excessive Sternheimer residuals: {residuals}')
        if name.startswith('epr'):
            ok, messages = check_output(kd, cv)
            if not ok:
                raise AssertionError('\n'.join(messages))
            values = parse_gtensor(kd.read_text())['delta_g_total']
            if serial is None:
                serial = values
            elif max(abs(a-b) for a,b in zip(values,serial)) > 0.02:
                raise AssertionError(f'Serial/MPI mismatch: {serial} vs {values}')
        else:
            # Compare using the repository's NMR checker and its physical tolerances.
            ref = dest / 'reference'
            ref.mkdir()
            shutil.copyfile(cv, ref / 'kdotp.out')
            subprocess.run([sys.executable, str(ROOT / 'check_nmr.py'), '--files',
                            str(kd), '--refdir', str(ref)], check=True)
        print(f'PASS {name}: kdotp vs covariant', flush=True)
    print('All kdotp integration checks passed.', flush=True)


if __name__ == '__main__':
    main()
