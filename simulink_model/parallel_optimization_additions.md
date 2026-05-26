# Parallel Optimization Additions

This document explains the parallel optimization changes added on top of the row-centered shading workflow. The goal is to make GA and MOGA runs work in either serial or parallel mode across the four main configurations:

- legacy fixed-axis
- legacy single-axis
- row-centered fixed-axis
- row-centered single-axis

## Why This Was Added

Large GA/MOGA runs spend most of their time evaluating many candidate designs through the Simulink model. MATLAB's Global Optimization Toolbox can distribute those candidate evaluations across MATLAB workers when `UseParallel` is enabled.

The main issue is that Simulink workers do not automatically inherit every initialized bus object and workspace variable from the client MATLAB session. If a worker tries to run the model with stale or missing bus definitions, the optimization can fail even when the same design works in serial.

## Main Implementation Changes

`agrivoltaic_wrapper.m` now builds the required Simulink bus objects from the actual `agriParams` and `agriVar` being simulated. This is important because the bus layout changes between legacy and row-centered mode, and between fixed-axis and single-axis runs.

The wrapper stores a lightweight bus signature in the worker base workspace. If a later simulation uses the same struct layout, the existing bus objects are reused. If the layout changes, the bus objects are rebuilt. This avoids relying on rerunning `agrivoltaics_variable_definition` inside MATLAB workers.

The optimizer scripts now use a `USE_PARALLEL_PROCESSING` toggle. When the toggle is `true`, the script starts a parallel pool only if one does not already exist. Existing user pools are not deleted automatically.

## Helper Files Added

- `functions/setup_parallel_pool.m`: starts or reuses a parallel pool and makes sure workers are in the model folder with the right MATLAB path.
- `functions/configure_agri_mode.m`: builds mode-safe `agriParams`, `agriVar`, lower bounds, and upper bounds for legacy/row-centered and fixed/single-axis runs.
- `functions/build_initial_population.m`: creates repeatable initial populations that include the current nominal design.
- `functions/build_tracking_slew_constraints.m`: creates the linear constraints used to enforce the tracking slew-rate limit.
- `functions/validate_tracking_slew.m`: checks whether optimized single-axis tracking angles satisfy the slew-rate limit.
- `parallel_validation_benchmark.m`: runs serial versus parallel validation and timing comparisons.

## Validation Workflow

Run this from the `simulink_model` folder:

```matlab
parallel_validation_benchmark
```

The benchmark writes a timestamped folder under:

```text
results/parallel_validation_*
```

Expected outputs include:

- `validation_summary.csv`
- `validation_report.md`
- `validation_raw.mat`
- `parallel_speedup.png`
- `objective_comparison.png`

The validation checks all four modes, compares serial and parallel wrapper outputs, confirms row-centered panel count stays at 550, and checks single-axis slew constraints.

## Practical Notes

Parallel mode is a correctness and scalability feature, not a guaranteed speedup for tiny tests. For short validation runs, worker startup and Simulink initialization overhead can make parallel slower than serial. Parallel should be most useful for larger GA/MOGA runs where the cost of many model evaluations is high enough to amortize the worker overhead.

Generated validation outputs should not be committed unless the group explicitly wants to archive a reference benchmark.
