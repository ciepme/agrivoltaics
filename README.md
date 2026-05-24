# agrivoltaics
ME599-001-bESSt_optimization_model_for_agrivotaic_systems

## Purpose
For University of Michigan's Mechanical Engineering 599 - 001 Course: Multidisciplinary Optimization.

## Scope
Berry production in Ann Arbor, MI using monocrystaline, monofacial silicon panels.

## Tools
Data from python's pvpy library and analysis using optimization tools through MATLAB/Simulink package

## Simulink Architecture

### Folder Structure

* ancillary_functions = optional programs that augment functionality such as SQP
* helper_functions = helper functions (e.g. wrappers and .m files called from mains)
* getNonDominated = 3rd party function to determine non-dominated set on a pareto front
* graphs = graphical outputs from other files
* modules = sections from main agrvoltaic model stored as library blocks
* parameter_data = idk
* results = data from run GA's and other model runs

### Core Files

### Ancillary Functions