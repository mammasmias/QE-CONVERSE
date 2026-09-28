"""Controlled CeO2 NMR experiment; never regenerates regression references.

A: original U path; B: promoted U with real occupations; C: complex occupations.
Builds use pinned sources. Run on Linux with QE 7.5, MPI and a Fortran compiler.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import urllib.request
import zipfile

BASE = '3bafb3f740039869bfbf9311f7ccbb82a5835e0f'
PATCH = '15b5f1689421249feb46f3970facbaf4bce44f07'
GIPAW = '3562150651c189a678b6bc8c9873bb62809189af'
PSEUDOS = {
    'Ce.pbe-tm-semi-dc-dipaw.UPF': '69fdcb02adc1f30d2b027105ca7d69215b3db1444afa0f4f89a03252b5e30c4c',
    'O.pbe-tm-gipaw-dc.UPF': '8c7c214e08ab43519522ab30061e06586868609b734f6b939a8982424fe9371c',
}
FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def build(root, variant, qe):
    sha = BASE if variant == 'A' else PATCH
    archive = zipfile.ZipFile(io.BytesIO(fetch(f'https://codeload.github.com/mammasmias/QE-CONVERSE/zip/{sha}')))
    archive.extractall(root)
    source = root / f'QE-CONVERSE-{sha}'
    if variant == 'B':
        path = source / 'src/new_nsg.f90'
        text = path.read_text()
        old = 'wg(ibnd,ik) * proj%k(off1,ibnd) * &\n                             CONJG(proj%k(off2,ibnd) * phase)'
        new = 'DBLE( wg(ibnd,ik) * proj%k(off1,ibnd) * &\n                             CONJG(proj%k(off2,ibnd) * phase) )'
        assert text.count(old) == 1, 'Diagnostic patch must match exactly once'
        path.write_text(text.replace(old, new))
        (root / 'B-only-change.txt').write_text(old + '\nREPLACED BY\n' + new)
    subprocess.run(['bash', './configure', f'--with-qe-source={qe}'], cwd=source, check=True)
    subprocess.run(['make', '-j2'], cwd=source, check=True)
    return source / 'bin/qe-converse.x'


def execute(command, text, directory, stem):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f'{stem}.in').write_text(text)
    with (directory / f'{stem}.out').open('w') as output:
        subprocess.run(command, input=text, text=True, cwd=directory, stdout=output,
                       stderr=subprocess.STDOUT, check=True, timeout=5400)
    result = (directory / f'{stem}.out').read_text()
    if 'convergence has been achieved' not in result or 'JOB DONE.' not in result:
        raise RuntimeError(f'Incomplete/nonconverged calculation: {directory}/{stem}.out')
    return result


def values(label, text):
    matches = re.findall(label + r'\s*(' + FLOAT + r'(?:\s+' + FLOAT + r')*)', text)
    if not matches:
        raise ValueError(f'Missing output field {label}')
    return [float(x.replace('D', 'E').replace('d', 'e')) for x in matches[-1].split()]


def pw_input(pseudo, use_u):
    text = f"""&control
 calculation='scf', prefix='ceo2', pseudo_dir='{pseudo}', outdir='./scratch/'
/
&system
 ibrav=2, a=5.409719944, nat=3, ntyp=2, ecutwfc=120,
 nosym=.true., noinv=.true.
/
&electrons
 conv_thr=1d-11, diagonalization='david', mixing_beta=0.2, electron_maxstep=200
/
ATOMIC_SPECIES
Ce 140.116 Ce.pbe-tm-semi-dc-dipaw.UPF
O 15.999 O.pbe-tm-gipaw-dc.UPF
K_POINTS automatic
4 4 4 0 0 0
ATOMIC_POSITIONS alat
Ce 0.00 0.00 0.00
O 0.25 0.25 0.25
O 0.75 0.75 0.75
"""
    if use_u:
        text += '\nHUBBARD (atomic)\nU Ce-4f 5.0\nU O-2p 1e-4\n'
    return text


def conv_input(atom, moment):
    return f"""&input_qeconverse
 prefix='ceo2', outdir='./scratch/', diagonalization='david',
 q_gipaw=0.01, dudk_method='covariant', mixing_beta=0.2,
 conv_threshold=1d-11, tr2=1d-11, m_0(1)={moment}, m_0_atom={atom},
 lhub_magnetization=.true.
/
"""


def experiment(args):
    root = args.root.resolve() / args.variant
    root.mkdir(parents=True, exist_ok=False)
    executable = build(root, args.variant, args.qe.resolve())
    pseudo = root / 'pseudo'
    pseudo.mkdir()
    for name, expected in PSEUDOS.items():
        data = fetch(f'https://raw.githubusercontent.com/dceresoli/qe-gipaw/{GIPAW}/pseudo/{name}')
        assert hashlib.sha256(data).hexdigest() == expected, name
        (pseudo / name).write_bytes(data)
    result = dict(variant=args.variant, source=BASE if args.variant == 'A' else PATCH,
                  gipaw_fixture=GIPAW, pseudos=PSEUDOS, cases=[], pool_checks=[])
    # Identical rank/pool layout and independently converged seeds for all variants.
    mpi = ['mpirun', '--oversubscribe', '-np', '4']
    for mode in ('U0', 'U5'):
        seed = root / mode / 'seed'
        execute(mpi + [str(args.qe.resolve() / 'bin/pw.x'), '-nk', '4'],
                pw_input(pseudo.as_posix(), mode == 'U5'), seed, 'pw')
        for atom in (1, 2):
            for amplitude in (1.0, 0.5):
                for sign in (1, -1):
                    case = root / mode / f'atom{atom}-m{sign * amplitude:+g}'
                    case.mkdir()
                    shutil.copytree(seed / 'scratch', case / 'scratch')
                    text = execute(mpi + [str(executable), '-nk', '4'],
                                   conv_input(atom, sign * amplitude), case, 'converse')
                    entry = dict(mode=mode, atom=atom, amplitude=amplitude, sign=sign,
                                 energy=values(r'!\s+total energy\s*=', text)[0],
                                 shift=values(r'Chemical shift\s*\(ppm\):', text),
                                 core=values(r'Core shift\s*\(ppm\):', text)[0])
                    for name, label in [('M', r'M_tot\s*='),
                                        ('M_hubbard', r'Delta_M_hubbard\s*='),
                                        ('imag', r'Complex Hubbard: max \|Im n\|\s*='),
                                        ('hermiticity', r'Complex Hubbard: onsite Hermiticity error\s*=')]:
                        entry[name] = values(label, text) if re.search(label, text) else None
                    if args.variant != 'A' and mode == 'U5':
                        assert 'on-site U initialized in the U+V representation' in text
                        assert entry['hermiticity'][0] < 1e-10
                        if args.variant == 'B':
                            assert entry['imag'][0] < 1e-12
                    result['cases'].append(entry)
                    (root / 'results.json').write_text(json.dumps(result, indent=2))
                    print(json.dumps(entry), flush=True)
                    # Only remove generated scratch after successful extraction; inputs/logs stay.
                    shutil.rmtree(case / 'scratch')
        if mode == 'U5':
            baseline = next(r for r in result['cases'] if
                            (r['mode'], r['atom'], r['amplitude'], r['sign']) == ('U5', 1, 1.0, 1))
            for pools in (1, 2):
                case = root / mode / f'atom1-pool{pools}'
                case.mkdir()
                shutil.copytree(seed / 'scratch', case / 'scratch')
                command = [str(executable)] if pools == 1 else ['mpirun', '--oversubscribe', '-np', '2', str(executable), '-nk', '2']
                text = execute(command, conv_input(1, 1.0), case, 'converse')
                shift = values(r'Chemical shift\s*\(ppm\):', text)
                result['pool_checks'].append(dict(pools=pools, shift=shift,
                    max_difference_ppm=max(abs(x-y) for x,y in zip(shift,baseline['shift']))))
                (root / 'results.json').write_text(json.dumps(result, indent=2))
                shutil.rmtree(case / 'scratch')
        shutil.rmtree(seed / 'scratch')


def compare(args):
    records = {}
    for file in args.root.rglob('results.json'):
        data = json.loads(file.read_text())
        assert data['variant'] not in records, 'Duplicate variant results'
        assert len(data['cases']) == 16, 'Incomplete experiment'
        records[data['variant']] = data
    assert set(records) == {'A', 'B', 'C'}, 'All three variants are required'
    central = {}
    lines = ['# CeO2: controlled Hubbard NMR comparison', '',
             'Printed CONVERSE valence response, ppm; no macroscopic/core convention matching to QE-GIPAW is assumed.',
             'The printed shift divides by |m|. Thus the odd response is (shift(+m)-shift(-m))/2.', '',
             '| Variant | U | Atom | |m| | Odd xx (ppm) | Even xx (ppm) |',
             '|---|---|---|---|---|---|']
    problems = []
    for variant, data in sorted(records.items()):
        assert len(data['pool_checks']) == 2, 'Missing serial/pool=2 checks'
        for check in data['pool_checks']:
            if check['max_difference_ppm'] > 2:
                problems.append(f"{variant}: pool={check['pools']} differs from pool=4 by {check['max_difference_ppm']:.4f} ppm")
        for mode in ('U0', 'U5'):
            for atom in (1, 2):
                for amplitude in (1.0, 0.5):
                    pair = {r['sign']: r for r in data['cases'] if (r['mode'], r['atom'], r['amplitude']) == (mode, atom, amplitude)}
                    assert set(pair) == {-1, 1}
                    odd = [(p-n)/2 for p,n in zip(pair[1]['shift'], pair[-1]['shift'])]
                    even = [(p+n)/2 for p,n in zip(pair[1]['shift'], pair[-1]['shift'])]
                    central[variant, mode, atom, amplitude] = odd
                    lines.append(f'| {variant} | {mode} | {atom} | {amplitude} | {odd[0]:.4f} | {even[0]:.4f} |')
                    if max(map(abs, even)) > 2:
                        problems.append(f'{variant}/{mode}/atom{atom}/m{amplitude}: even residual > 2 ppm')
    lines += ['', '## Differences and diagnostic gates', '']
    for mode in ('U0', 'U5'):
        for atom in (1, 2):
            for amplitude in (1.0, 0.5):
                a, b, c = [central[v, mode, atom, amplitude] for v in 'ABC']
                ab = max(abs(x-y) for x,y in zip(a,b))
                bc = max(abs(x-y) for x,y in zip(b,c))
                lines.append(f'- {mode}, atom {atom}, |m|={amplitude}: max|B-A|={ab:.4f} ppm; max|C-B|={bc:.4f} ppm.')
                if ab > 2:
                    problems.append(f'{mode}/atom{atom}/m{amplitude}: A-B inequivalence > 2 ppm')
                if mode == 'U0' and bc > 2:
                    problems.append(f'U0/atom{atom}/m{amplitude}: B-C differs > 2 ppm')
    for variant in 'ABC':
        for mode in ('U0', 'U5'):
            for atom in (1, 2):
                delta = max(abs(x-y) for x,y in zip(central[variant,mode,atom,1.0],central[variant,mode,atom,0.5]))
                if delta > 2:
                    problems.append(f'{variant}/{mode}/atom{atom}: amplitude dependence {delta:.4f} ppm > 2 ppm')
    lines += ['', 'A-B must agree before attributing B-C to imaginary occupations. C-B at U5 is measured, not required to vanish.',
              'This is a fixed-grid diagnostic, not a cutoff/k-point convergence study or proof of the complete magnetic functional.', '',
              '## Outcome', ''] + (problems or ['All diagnostic gates passed.'])
    report = '\n'.join(lines) + '\n'
    (args.root / 'comparison.md').write_text(report)
    print(report)
    summary = os.environ.get('GITHUB_STEP_SUMMARY')
    if summary:
        with open(summary, 'a') as f:
            f.write(report)
    if problems:
        raise SystemExit('Diagnostic differences require investigation; references were not changed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run', 'compare'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--variant', choices=list('ABC'))
    parser.add_argument('--qe', type=Path)
    args = parser.parse_args()
    if args.action == 'run':
        if not args.variant or not args.qe:
            parser.error('run requires --variant and --qe')
        experiment(args)
    else:
        compare(args)
