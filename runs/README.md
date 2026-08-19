# Run artifacts

Per-condition validation logs and Klempt A-matrix environment configs emitted by
the growth/biofilm runs (`p23_*`, `p30_*`, `p31_*`). These are **evidence
outputs**, not inputs — nothing in the codebase reads them by path (verified).
Kept for provenance; regenerate from the corresponding run scripts.

## Computational cost

Abaqus writes the cost evidence the thesis needs into plain-text job files that
no ODB tooling is required to read:

| File | Carries |
|---|---|
| `<job>.sta` | increment history — increments, attempts, cutbacks, equilibrium and severe-discontinuity iterations |
| `<job>.msg` | `JOB TIME SUMMARY` — wallclock and CPU seconds; warning/error counts |
| `<job>.dat` | model size — elements, nodes, DOF |

Extract them with [`extract_abaqus_cost.py`](../extract_abaqus_cost.py) (stdlib
only; no Abaqus licence, no ODB):

```bash
# scan a scratch directory and commit the result here
python3 extract_abaqus_cost.py /path/to/abaqus/scratch \
        -o runs/abaqus_cost.json --markdown
```

`--markdown` prints a table ready to paste into the cost section. Derived
columns: iterations per increment, seconds per equilibrium iteration, and
microseconds per element-iteration — the last being the figure to quote when
comparing UMAT cost against a reference material.

**A value that could not be read is `null`, never `0`.** A job with only a
`.sta` yields iteration counts with `wallclock_s: null`, and the script says on
stderr which jobs had no timing. Do not fill those in by hand.
