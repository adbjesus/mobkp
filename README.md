# mobkp

### Introduction

This repository contains several (anytime) algorithms to solve the
Multi-Objective Binary Knapsack Problem (MOBKP). It is made up of a solver
`mobkp`, which can be found on the `apps/` directory, and an header only
library under the namespace `mobkp`, which can be found on the `include/`
directory.

### Dependencies

This library depends on

- [mooutils](https://github.com/adbjesus/mooutils) for multi-objective
  optimization utilities, such as quality indicators and solution sets.
- [apm](https://github.com/adbjesus/apm) for a theoretical anytime performance
  model that guides one of the algorithms.
- [glpk](https://www.gnu.org/software/glpk)
- [fmt](https://github.com/fmtlib/fmt)
- [CL11](https://github.com/CLIUtils/CLI11)
- [Boost](https://boost.org)

### Compiling/Installing

You can use `cmake` to compile/install the library and/or solver.
