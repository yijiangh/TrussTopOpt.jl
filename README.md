# TrussTopOpt.jl

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/yijiangh/TrussTopOpt.jl.svg?branch=master)](https://travis-ci.org/yijiangh/TrussTopOpt.jl)
[![codecov](https://codecov.io/gh/yijiangh/TrussTopOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yijiangh/TrussTopOpt.jl)

Truss topology optimization experiments using [TopOpt.jl](https://github.com/mohamed82008/TopOpt.jl) and [JuAFEM.jl](https://github.com/KristofferC/JuAFEM.jl).

## Documentation

<!-- [![Documentation](https://img.shields.io/badge/doc-latest-blue.svg)](https://yijiangh.github.io/TopOpt.jl/dev) -->
[Dev doc](https://docs.google.com/document/d/1SbHCct2B0aTH4drSAJX5HGYlxkyCiAznsFe5GMzzCyw/edit?usp=sharing)

## Installation

In Julia v1.0+ you can install packages from the Pkg REPL (press `]` in the Julia
REPL to enter `pkg>` mode):

```
pkg> add https://github.com/yijiangh/JuAFEM.jl.git
pkg> add https://github.com/yijiangh/Tensors.jl.git
pkg> add https://github.com/mohamed82008/VTKDataTypes.jl#master
pkg> add https://github.com/mohamed82008/KissThreading.jl#master
pkg> add https://github.com/yijiangh/TopOpt.jl#master
```

which will track the `master` branch of the package.

To load the package, use

```julia
using TrussTopOpt
```