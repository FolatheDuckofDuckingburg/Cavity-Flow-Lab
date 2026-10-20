# Validation protocol (expanded)

See [VALIDATION.md](VALIDATION.md) for the main project validation plan. This document describes the separate gates and their limits.

| Gate | Evidence required | Does not prove |
|---|---|---|
| Software invariants | Passing tests for shapes, configuration, finite fields, and step behavior | Physical accuracy |
| Projection | Independently checked discrete operators, pre/post divergence, pressure residual | Benchmark accuracy |
| Reference comparison | Verified Ghia data, sampling method, RMSE and maximum error | General accuracy at all grids/Re |
| Grid/time sensitivity | Comparable physical time or justified steady state, multiple resolutions and time steps | Formal order without appropriate analysis |
| Performance | Repeated timings on same hardware/configuration with correctness checks | Better physics merely because it is faster |

## Reporting checklist

For every result record commit SHA, software environment, grid, Reynolds number, time-step policy, steps, simulated time, pressure iteration count, residual, pre/post-projection divergence, centreline error metrics where reference data are available, runtime distribution, and limitations.

The current scripts provide an initial reproducibility layer, not every quantity in the checklist. The projection operator should be independently audited, and reference data must be verified before automated benchmark acceptance tests are added.
