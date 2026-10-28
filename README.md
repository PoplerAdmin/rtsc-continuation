# rtsc continuation work

Continuation of the `rtsc` superconductivity-screening project from the 2026-10-27 handoff.

| Path | Contents |
|---|---|
| `HANDOFF.md`, `STARTER_PROMPT.txt`, `rtsc_handoff_2026-10-27/` (+ `.tar.gz`) | Handoff bundle as received |
| `verify_results_2026-10-27/` | Step-1 verification of the bundle (Lean and Python certificate checks, demo run) |
| `reports/` | Benchmark, convergence-protocol and tetrahedron reports; `reports/data/` holds their records |
| `bench_tools/` | Analysis, EPW driver and the `ephproto` parsers/admission code, with tests |
| `run_scripts/scratchpad_2026-10-27/` | The working scripts that produced the results, kept verbatim apart from the path changes below |

## Anonymized paths

Absolute paths were rewritten before publication. Nothing else was changed.

- The WSL home directory of the run machine appears as `/home/user`, and the Windows profile as `C:\Users\user`.
- References to the agent's temporary scratch folder point to `run_scripts/scratchpad_2026-10-27`, where copies of those scripts now live. Leftover temporary-folder session IDs appear as `<session>`.

To rerun a script, replace `/home/user` with your own paths. Scripts that embed paths in single quotes, Python or QE inputs need the literal path; `$HOME` will not expand there.

`config_hash` values in the run records were computed on the original paths. `rtsc_run_config.toml` includes the DFT binary and pseudopotential directories, so recomputing the hash from the anonymized file gives a different value. That hash only identifies a run configuration on one machine. Certificates do not depend on it: `rtsc recheck` re-derives every verdict from the α²F spectrum and the certificate parameters, and those were not changed.
