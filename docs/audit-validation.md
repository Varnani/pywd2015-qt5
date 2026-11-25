# PyWD2015 Qt5 audit corrections

Baseline: `756aab0d6af826ddd857fb850ae226c6ff0d0bda`. Work is on local `fix/audit-correctness`; no remote publication or new GUI executable.

## Implemented findings

| Finding | Correction and evidence |
| --- | --- |
| B01 | Plaintext/LaTeX export reads output column 4; fractional Kelvin and sigma preserved. `test_export.py`. |
| B02 | Child nonzero exit, missing output and malformed numerical tables fail before publication; previous shared output remains untouched. `test_process.py`, `test_lifecycle.py`. |
| B03 | Worker catches exceptions, prepares runtime before notification, retains cancellation before spawn, and UI cleanup waits for QThread.finished. Repeated abort/restart and error restoration tests. |
| B04 | Each IO owns a fresh temporary job directory; executable path is absolute, auxiliary `.dat`/`.tab` files copied, stdout/stderr remain associated with the job. Separate application processes preserve distinct inputs. |
| B05 | Blocking child work runs on QThread while a modal progress dialog pumps the Qt event loop. Timer, cancellation and timeout tests; no plotting/widget changes run in that thread. |
| B06 | Only fractional mantissa zeroes are trimmed; integer seed digits and exponent digits survive; F fields use fixed decimal syntax. `test_numerics.py`. |
| B07 | Per-observation cycle bounds with exact coordinate filtering; large-JD endpoint rounding handled. Brute-force reference tests. |
| B08 | Logical DC equations derive width from final fitted parameter count, jdphs and ifoc. Fixed Fortran fields parsed across any number of physical lines. Equations, stats and dimensions belong to the final solution; partial final tables cannot borrow earlier data. Synthetic 38/41 and 20/23 boundaries plus real native fixtures. |
| B09 | Shared two-dimensional observation reader checks shape, finite data, selected columns and nonnegative weights; two-column eclipse data gets unit weights. Failed loads retain old state. |
| B10 | Historical Julian/Gregorian policy in both directions, including century/leap and transition tests; invalid input remains inside the GUI error boundary. No timezone/time-scale conversion. |
| B11 | Explicit shared phase control remains enabled in both position/Roche modes; period is no longer used as phase. `.ui` sources and generated Python updated. |
| B12 | Independent project validation checks all 165 existing static UI fields, dynamic curves/spots and observation paths before applying; signal blockers protect bulk application. KEEP selections restored after constraints. Atomic save. |
| B13 | Project, spot and application INI use raw text interpolation=None. UTF-8 new saves; UTF-8 then original machine-locale fallback for legacy files. Literal percent/Turkish path tests. |
| B14 | Sampling bounds/step checked before division and LC input serialization; empty/malformed motor frames rejected. |
| M01 | CPython 3.12 source support documented; old Python/sudo-pip instructions removed; correct PyQt5 name and venv setup. |
| M02 | Hash-locked runtime and development dependencies; Windows/Ubuntu offscreen CI configured; real WD fixtures with provenance and golden parser tests. |
| M03 | `tools/release_manifest.py` records source commit/dirty state/source hashes, engine+auxiliary hashes, dependency versions and SBOM hash. No replacement application binary released. |
| M04 | shell=False argument list with absolute executable, hidden Windows child; real binaries tested from spaced/Turkish executable directories. |
| M05 | Empty history does not open export dialog or index its first row. |
| M06 | Resource paths derive from module __file__; settings use the user's Qt config directory, with legacy-file fallback. |

## Native engine reference

Official v1.0.0 `wd-precompiled.zip` SHA-256: `5bfa347dfcaf858e5cf3c3bcd8bbc9adf2ce553906df5a32e0bd519f49279116`.

Exact engine/source/data and fixture hashes are in `tests/fixtures/wd2015-win64/provenance.json`. Binaries remain in ignored `.local/engines`; source delivery includes outputs and provenance, not redistributed engine executables.

Real Windows LC: 21 velocity observations of an eccentric detached system, e=0.2, inclination=85 degrees, period=2 days, a=5 solar radii, q=0.8, temperatures 6500/6000 K, potentials 5. Real DC fits semimajor axis and systemic velocity to those observations with +0.5 km/s offset. One-iteration output: a=5.0000011602, gamma=0.4999956211 km/s. Two iterations: gamma=0.5000006220 km/s. Real fixtures cover phase/time and ifoc=0/1; unsupported executable variants are not inferred to share this contract.

Regenerate optional integration checks on Windows:

```powershell
$env:PYWD_ENGINE_DIR = 'absolute path to extracted wd-precompiled'
.venv/Scripts/python.exe -m pytest -q
```

## Validation limits

- Windows CPython 3.12.11 and pinned dependency set tested offscreen. Linux CI is supplied but was not remotely executed; no macOS support claim.
- Native golden outputs establish reproducibility of this executable/data set and parser, not independent scientific certification of the WD physics. Large 38/41 and upstream 20/23 wrap boundary cases are synthetic; real native fixtures use two fitted parameters. Additional motor builds and large real physical solutions need their own fixtures.
- No interactive user desktop was driven. Physical close/resize/plot acceptance remains a user-facing follow-up; automated Qt cancellation/restart and timer tests run offscreen.
- Python advisory scan covers the locked Windows runtime package set. It reports no known advisories in that set; this does not certify native Qt, Fortran code, or the old packaged application EXE.
- Existing invalid-escape SyntaxWarnings in original LaTeX strings are outside this correction scope. No PyQt6 migration or physics formula rewrite.
- Job files are temporary and only the most recent IO is retained for the viewer; export durable results before exiting. Logs survive for the latest failed job while the process is open.

## Verification results

Final Windows run: `PYWD_ENGINE_DIR=<verified official engine directory> python -m pytest -q --tb=short` — **105 passed, 0 skipped**, 45.90 seconds. The environment was synchronized with `uv pip sync --require-hashes requirements-dev.lock`; 46 locked development/runtime packages matched. This includes real LC/DC execution from Turkish/spaced executable directories.

`python -m compileall -q src gui tests tools` and `git diff --check` passed. Runtime `pip-audit -r requirements.lock --no-deps --disable-pip` returned no known vulnerabilities; JSON and CycloneDX output are included. An independent review found and verified fixes for invalid-date Qt exceptions, KEEP restoration, legacy encodings, curve/spot ranges, final-iteration table scope and dialog cleanup. CI workflow is configured but not run remotely.
