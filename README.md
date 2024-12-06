# Phaseless auxiliary-field quantum Monte Carlo method for cavity-QED matter systems
[![DOI](https://zenodo.org/badge/899588988.svg)](https://doi.org/10.5281/zenodo.14289884)

This is the data repository for our paper [Phaseless auxiliary-field quantum Monte Carlo method for cavity-QED matter systems](https://doi.org/10.48550/arXiv.2410.18838).

The data is in JSON format in the `data` directory.

To regenerate the plots, you need Julia. Running

```bash
cd scripts
julia --project -e "using Pkg; Pkg.instantiate()"
julia --project make_plots.jl
```

should install all further dependencies and write plots to the `figs` directory.
