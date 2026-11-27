# Audit correction plan

Source: user-provided PyWD2015 Qt5 audit, baseline 756aab0d6af826ddd857fb850ae226c6ff0d0bda.

- [x] B01/B06: reproduce export and fixed-width serialization errors, preserve units and precision.
- [x] B02–B05/M04: isolate job files, validate child exit/output, cancellation/timeout, safe Qt worker lifecycle and responsive LC execution.
- [x] B07/B10: per-observation alias cycles and consistent historical calendar conversion.
- [x] B08: solution-scoped logical DC equation records, verify native engine contract.
- [x] B09/B12/B13: validated observation IO, staged project model, raw INI and atomic save.
- [x] B11/B14/M05: usable Roche phase control, checked sample intervals and empty history.
- [x] M01–M03/M06: pinned dependencies, Windows/Linux CI, resource/config locations, release provenance tooling.
- [x] Final: run full regression suite, dependency audit, native LC/DC fixture checks, independent review; report any unverified acceptance criteria.

Keep the MIT license and existing physics/formats. The 2026-11-28 follow-up authorizes pushing the correction branch and opening an upstream PR; release and desktop interaction are outside this request. Each nontrivial correction gets a failing regression test before implementation. Native engine validation must be distinguished from controlled process fixtures.

Ruling: Roche uses an explicitly labeled selectable phase, separate from orbital period; reuse a shared control outside either mode group.

Validation boundaries and externally unexecuted checks are recorded in audit-validation.md; CI/native scientific certification are not claimed complete.

- [x] PR follow-up: reproduce/fix final-iteration curve statistics, independently review lifecycle, and generate real wide native phase/time fixtures.
- [x] Final expanded Windows suite: 112 passed with native engines; repeat locally and inspect Windows/Linux CI during PR submission.
