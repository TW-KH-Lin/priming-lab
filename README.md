# PrimingLab

PrimingLab is a browser-based implementation of a sequential vesicle priming and short-term plasticity model. It runs entirely on the visitor's device: no Python installation, server-side computation, account, or data upload is required.

## First release

- Four priming-model branches (ES ⇄ LS ⇄ TS, optional TSL and ERS)
- PNAS 2022 Figure 4 and JP286282 example presets
- Editable release, priming, calcium, facilitation, and depression parameters
- Multi-frequency simulation with fixed-step RK4 integration
- Normalized release, pool occupancy, and release-probability plots
- CSV export
- Responsive desktop and mobile interface

## Run locally

Serve the repository with any static web server, for example:

```bash
python3 -m http.server 4173
```

Then open http://127.0.0.1:4173.

## Scientific scope

The web simulator mirrors the homogeneous fixed-step equations in the accompanying Python model. It is intended for research exploration, teaching, and qualitative comparison. Confirm quantitative results against the reference Python implementation before publication.

Primary reference: [PNAS 2022, DOI 10.1073/pnas.2207987119](https://doi.org/10.1073/pnas.2207987119).

## Privacy

Simulation and CSV generation happen locally in the browser. PrimingLab does not upload experimental data.

## License

No software license has been selected yet. The repository owner should choose a license before inviting redistribution or modification.
