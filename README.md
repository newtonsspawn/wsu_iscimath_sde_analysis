# Stochastic Differential Equation Analyses

A research repository for **numerical stochastic differential equations (SDEs)**, focused on **strong vs. weak convergence**, **method comparison**, and **reproducible Monte Carlo experiments**.

The primary workflow lives in notebooks and modular research code, with a project layout that supports gradually moving validated components into a more “production” style structure.

---

## Overview

This repo is currently a lab studying convergence for SDE solvers. It is designed to make it easy to:

- Implement canonical discretization schemes (e.g., Euler–Maruyama, Milstein, and Stratonovich-oriented methods).
- Benchmark them on SDEs with known properties/closed forms (starting with GBM as a baseline).
- Measure **strong error**, **weak error**, and **runtime vs. accuracy** tradeoffs.
- Investigate practical pitfalls in Monte Carlo estimation (e.g., sampling noise dominating weak error curves).
- Explore how modeling assumptions (notably the increment/noise model) affect observed convergence behavior.

---

## Goals

### Core goals
- **Convergence validation**: Empirically verify expected **strong** and **weak** convergence behavior for standard schemes.
- **Method comparison**: Compare solvers by accuracy, stability, and runtime scaling.
- **Fair Monte Carlo diagnostics**: Use best practices like common random numbers / nested increments when appropriate.
- **Reproducibility**: Ensure experiments can be rerun consistently via a pinned environment.

### Secondary goals
- Expand beyond baseline examples to additional SDE families (multi-dimensional systems, non-Lipschitz diffusions, etc.).
- Promote reusable components from `/research` into `/src` once stabilized.
- Add tests that validate convergence and numerical properties (not just unit-level correctness).

---

## What’s inside

### Notebooks (research / experiments)
You’ll find most analysis work under:

- `research/analyses/notebooks/`

Current analyses include:
- **Weak vs. strong convergence experiments** (Euler–Maruyama vs. Milstein)
- **Itô vs. Stratonovich** demonstrations and drift conversion intuition

> Notebooks are the focal point of the repo right now.

### Modular research code
- `research/analyses/modules/` — reusable building blocks for experiments (e.g., stepping functions, error calculators, plotting utilities).

### Papers / reading notes
- `research/papers/` — papers and supporting material (as the folder grows).

### Production area (placeholder for stabilized code)
- `src/` — intended for production-ready implementations once research code is validated and refactored.
- `tests/` — intended for automated testing (including eventual numerical regression tests).

### Repository Structure

```
/wsu_iscimath_sde_analysis
├── /research               # Research environment
│   ├── /analyses           # Analyses not ready for production
│   │   ├── /modules        # Modular components for SDE analyses
│   │   ├── /notebooks      # Scratch Jupyter notebooks for SDE analyses
│   │   └── main.py         # Main entry point for modular analyses
│   ├── /papers             # Research papers
│
├── /src                    # Production environment
│   ├── /empty              #
│   │   ├── /empty          # 
│   │   └── null            #
│   
├── /tests                  # Test environment
│   ├── /empty              #
│   │   ├── /empty          # 
│   │   └── null            #
│
├── CHANGELOG.md            # Currently empty
├── CONTRIBUTING.md         # Currently empty
├── README.md
└── setup_instructions.md
```

---

## How to run

This repo uses **Pixi** to create a reproducible Python environment. Full setup steps live in: **[setup_instructions.md](setup_instructions.md)**

### Quick start:

1. Install Pixi: ```% pixi install```
2. Run a command inside the environment: ```% pixi run python -V```
3. Launch JupyterLab: ```% pixi run jupyter lab```
4. Run tests: ```% pixi run test```


---

## Roadmap

### Near-term
- [X] Add Pixi environment files + tasks (`pixi.toml`, `pixi.lock`) and wire into setup docs
- [ ] Consolidate common solver / error-measurement utilities in `research/analyses/modules/`
- [ ] Standardize experiment configuration (step sizes, Monte Carlo sample sizes, seeds)
- [ ] Add trajectory diagnostics for weak convergence (common random numbers / shared paths)
- [ ] Add automated numerical regression tests to validate the convergence pipeline and detect breakdowns
- [ ] Add “weak error with noise floor” plots to prevent misinterpreting Monte Carlo variance as convergence

### Mid-term
- [ ] Add more benchmark SDEs beyond GBM
- [ ] Add solver interfaces and shared experiment runners (so notebooks become thinner)
- [ ] Build numerical regression tests (e.g., “strong slope is near expected order” within tolerance)

### Longer-term
- [ ] Promote stabilized implementations into `/src`
- [ ] Add CI-friendly test suite for reproducibility and numerical sanity checks
- [ ] Produce a small set of “reference analyses” notebooks that serve as canonical examples

---

## References

Researcxh papers:
- [A. Rössler, Second order Runge–Kutta methods for Itô stochastic differential equations, SIAM J. Numer. Anal. 47 (2009), no. 3, 1713–1738.](research/papers/2nd Order Runge-Kutta Methods for ITO SDEs.pdf)

General references for the topics explored here:
- Kloeden, P. E., & Platen, E. *Numerical Solution of Stochastic Differential Equations.*
- Higham, D. J. (2001). *An Algorithmic Introduction to Numerical Simulation of Stochastic Differential Equations.*
- Øksendal, B. *Stochastic Differential Equations: An Introduction with Applications.*

---

## Credits

- Contributors: see [CONTRIBUTING.md](CONTRIBUTING.md)







