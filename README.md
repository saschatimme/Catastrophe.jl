# Catastrophes in Elastic Tensegrity Frameworks

This code accomplishes the article `Catastrophe in Elastic Tensegrity Frameworks`.

## Installation

```jula
using Pkg
Pkg.add("https://github.com/saschatimme/Catastrophe.jl.git")
```

## Usage

We provide two demos.

**Zeeman's catastrophe machine**
```julia
using Catastrophe
Zeeman.start_demo()
```

**Elastic four-bar framework**
```julia
using Catastrophe
FourBars.start_demo()
```
