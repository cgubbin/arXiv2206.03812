# Arxiv2206.03812

This repository contains code used to generate all figures in [Arxiv2206.03812](https://arxiv.org/abs/2206.03812). The simulation is carried out using [Julia](https://julialang.org/) which can be installed from source using the [distributed installer](https://julialang.org/downloads/) for your hardware platform.

When Julia is installed a full calculation, generating all figures from the manuscript for a single film thickness can be carried out running

```bash
julia --project=. scripts/FullRun.jl
```

from the top-level folder of the repository. This will download and pre-compile all dependancies before running the simulation.
