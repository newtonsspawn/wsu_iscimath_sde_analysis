# Pixi Environment Setup Instructions

This repo uses **Pixi** to create a reproducible Python environment. Instructions for setting up Pixi can be found [here](https://github.com/pixi-project/pixi).

---

## Environment Setup

1. Install Pixi -- Follow the official Pixi installation instructions for your OS, then confirm it works:

```bash
pixi --version
```

2. Create the environment -- From the `wsu_iscimath_sde_analysis/` directory:

```bash
pixi install
```

3. Run commands inside the environment -- Run a command in the Pixi environment like this:

```bash
pixi run python -V
```

## Available Tasks

Run Pixi tests:

```bash
pixi run test
```

```bash
pixi run pytest
```

Launch Jupyter Labs to run notebooks:

```bash
pixi run jupyter lab
``` 

#### Notes
  - Dependency definitions live in `pixi.toml`.
  - The exact resolved environment is locked in `pixi.lock` (commit this for reproducibility).
  - If you’re using PyCharm, set the project interpreter to the Python executable inside the Pixi environment for `wsu_iscimath_sde_analysis/`.