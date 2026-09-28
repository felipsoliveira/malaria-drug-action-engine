# Contributing to Malaria Drug Action Engine

Thank you for considering a contribution. MDAE is a scientific-computing project, so a change is only complete when it is both **software-correct** and **scientifically auditable**.

## Ways to contribute

Useful contributions include:

- compartment/exposure models for *Plasmodium falciparum*;
- experimentally grounded parameters with primary-source provenance;
- validation datasets and falsifiable benchmark cases;
- numerical-analysis improvements;
- tests for scientific invariants and edge cases;
- documentation, examples, and reproducibility fixes.

For substantial scientific changes, opening an issue first is encouraged so assumptions and validation criteria can be discussed before implementation.

## Development setup

MDAE requires Python 3.10+.

```bash
git clone https://github.com/felipsoliveira/malaria-drug-action-engine.git
cd malaria-drug-action-engine
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest -q
```

## Scientific contribution requirements

When a pull request changes model behaviour, parameters, or benchmark logic, please include:

1. **Scientific question** — what hypothesis or behaviour is being tested?
2. **Assumptions** — units, boundary conditions, approximations, and model scope.
3. **Provenance** — primary references for measured or literature-derived parameters.
4. **Uncertainty** — measured, derived, fitted, or assumed values should be distinguishable.
5. **Validation** — a reference result, analytical oracle, held-out data, or another explicit validation target.
6. **Failure criteria** — when practical, define PASS/FAIL or stopping criteria before inspecting the final result.
7. **Tests** — add or update automated tests for the behaviour being changed.

Negative results are valid results. Do not tune acceptance criteria after seeing an outcome simply to make a benchmark pass.

## Code-review checklist

Reviewers should check more than whether the code executes:

- Is the requested behaviour implemented?
- Are changes scoped and understandable?
- Do existing and new tests pass?
- Are units and numerical tolerances explicit?
- Are floating-point comparisons appropriate?
- Are physical or mathematical constraints preserved?
- Are parameter sources and assumptions documented?
- Is the workflow reproducible from a clean environment?
- Does the conclusion match the evidence without overclaiming?

## Pull requests

Keep pull requests focused. In the description, explain what changed, why it changed, how it was tested, and any scientific limitations.

The repository uses CI to run the test suite on supported Python versions. A failing CI run should be investigated rather than bypassed unless the failure is demonstrably unrelated and documented.

## Style

Prefer clear scientific code over clever code:

- explicit units and variable names;
- small, testable functions;
- deterministic behaviour where possible;
- fixed random seeds for benchmark workflows when stochasticity is not itself under study;
- comments for scientific assumptions, not obvious syntax.

## Reporting problems

For bugs or scientific inconsistencies, open an issue with a minimal reproducible example, expected behaviour, observed behaviour, environment information, and relevant references when applicable.
