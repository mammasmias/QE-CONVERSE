"""Run actual U/U+V binaries, checking imaginary occupations and pool agreement.

PW and QECONVERSE must name executables. Uses temporary run directories;
prints their location and keeps all inputs/outputs for inspection.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'


def numbers(pattern, text):
    matches = re.findall(pattern, text)
    if not matches:
        raise AssertionError(f'Missing output: {pattern}')
    fields = matches[-1] if isinstance(matches[-1], tuple) else (matches[-1],)
    return [float(x.replace('D', 'E').replace('d', 'e')) for x in fields]


def run(command, input_text, directory, stem):
    (directory / f'{stem}.in').write_text(input_text)
    with (directory / f'{stem}.out').open('w') as output:
        subprocess.run(command, input=input_text, text=True, cwd=directory,
                       stdout=output, stderr=subprocess.STDOUT, check=True, timeout=900)
    text = (directory / f'{stem}.out').read_text()
    assert 'convergence has been achieved' in text, f'SCF failed: {directory}/{stem}.out'
    assert 'JOB DONE.' in text, f'Incomplete run: {directory}/{stem}.out'
    return text


def main():
    here = Path(__file__).resolve().parent
    pw = shlex.split(os.environ['PW'])
    converse = shlex.split(os.environ['QECONVERSE'])
    mpi = shlex.split(os.environ.get('MPIEXEC', 'mpirun --oversubscribe'))
    root = Path(tempfile.mkdtemp(prefix='qe-complex-hubbard-'))
    print(f'Test inputs and outputs: {root}', flush=True)
    pw_input = (here / 'CO+' / 'pw_scf.in').read_text()
    pw_input = pw_input.replace("pseudo_dir = '../test_pseudos/'",
                                f"pseudo_dir = '{(here / 'test_pseudos').as_posix()}'")
    pw_input = pw_input.replace('    nbnd = 6', '    nbnd = 8\n    nosym = .true.\n    noinv = .true.')
    pw_input = pw_input.replace('1 1 1 0 0 0', '4 1 1 0 0 0')
    pw_input = pw_input.replace('conv_thr = 1e-8', 'conv_thr = 1e-10')
    conv_input = (here / 'CO+' / 'gtensor_1.in').read_text()
    conv_input = conv_input.replace('conv_threshold = 1d-8', 'conv_threshold = 1d-10')
    for mode in ('U', 'UV'):
        # C is site 2, O is site 1. Both have one p manifold.
        hubbard = '\nHUBBARD (ortho-atomic)\nU C-2p 2.0\nU O-2p 1.0\n'
        if mode == 'UV':
            hubbard += 'V C-2p O-2p 2 1 0.5\n'
        baseline = None
        for pools in (1, 2, 4):
            directory = root / f'{mode}-pool{pools}'
            directory.mkdir()
            run(pw, pw_input + hubbard, directory, 'pw')
            command = converse if pools == 1 else mpi + ['-np', str(pools)] + converse + ['-nk', str(pools)]
            text = run(command, conv_input, directory, 'converse')
            promoted = 'on-site U initialized in the U+V representation' in text
            assert promoted == (mode == 'U'), f'Wrong Hubbard dispatch in {directory}'
            imag = numbers(r'Complex Hubbard: max \|Im n\| =\s*(' + FLOAT + ')', text)[0]
            herm = numbers(r'Complex Hubbard: onsite Hermiticity error =\s*(' + FLOAT + ')', text)[0]
            assert imag > 1e-10, f'Imaginary occupations lost: {directory}'
            assert herm < 1e-10, f'Non-Hermitian occupations: {directory}'
            e = numbers(r'!\s+total energy\s*=\s*(' + FLOAT + ')', text)[0]
            m = numbers(r'M_tot\s*=\s*(' + FLOAT + r')\s+(' + FLOAT + r')\s+(' + FLOAT + ')', text)
            if baseline is None:
                baseline = (e, m)
            else:
                assert abs(e - baseline[0]) < 1e-6, f'Pool-dependent energy: {directory}'
                assert max(abs(a-b) for a, b in zip(m, baseline[1])) < 1e-4, f'Pool-dependent M: {directory}'
            print(f'PASS {mode} npool={pools}, max|Im n|={imag:.6g}', flush=True)


if __name__ == '__main__':
    main()
