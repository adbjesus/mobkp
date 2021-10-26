# Anytime-MOBKP

### Introduction

This repository contains several anytime algorithms to solve the
Multi-Objective Binary Knapsack Problem. It is made up of a solver
`mobkpsolve`, which can be found on the `apps/` directory, and an header
only library `mobkp`, which can be found on the `include/` directory.

If you use this code for your publications please cite it as: TODO

### Dependencies

The library depends on

- [mooutils](https://github.com/adbjesus/mooutils) for multi-objective
  optimization utilities, such as quality indicators and solution sets.

The solver depends on the library (and all its dependencies) and also
on:

- [fmt](https://github.com/fmtlib/fmt) for formatting and printing.
- [CLI11](https://github.com/CLIUtils/CLI11) for the cli interface.

### Compiling/Installing

You can use `cmake` to compile/install the library and/or solver.
