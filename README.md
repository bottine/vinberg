# 𝐕𝐢𝐧𝐧𝐲 𝐭𝐡𝐞 𝐜𝐨𝐱𝐛𝐨𝐲

### *Learning to herd Coxeter diagrams in the Hyperbolic Plains since 2021*

***

# Port of [aperep/vinberg-algorithm](https://github.com/aperep/vinberg-algorithm) to julia with some changes


**Warning.** The code is still being tested, and surely has many bugs!

## Description

* The finite volume check for Coxeter diagrams is implemented “from scratch” using Guglielmetti's thesis as theoretical resource and [`CoxIter`](https://github.com/rgugliel/CoxIter) as a comparison point (plus a majority of tests from there); thank you!
* The main bulk of Vinberg's algorithm (and associated procedures) is ported (modulo errors introduced in the process) from [B&P](https://github.com/aperep/vinberg-algorithm).
* The code uses static bitsets for hopefully efficient storage of small integer sets.
  The code for this comes from [here](https://raw.githubusercontent.com/chethega/StaticArrays.jl/fb0350012f01db4021d60906357e949333ec5d93/src/SBitSet.jl).
  **TODO.** Contact the author to ask if we can use their code, and under what licence (see [here](https://github.com/chethega/StaticBitsets.jl/issues/1)).

## How to use

* If not already installed, run `julia install_packages.jl` to install the packages used in the code.
* Go to the root folder of the project (having subfolder `lattices`,`graphs`,`src`).
* Launch julia with `julia -t auto` to enable multithreading (might be useful, but on small examples does't look like it is)
* call `include("src/vinberg.jl")`.
* calling `toggle(true)` or `toggle(false)` respectively enables/disables a majority of `asserts` in the code.

Now everything needed is loaded.

* To run the Vinberg Algorithm on a lattice given by a matrix `M`, call `Vinberg_Algorithm(M)`.
* Some matrices are accessible through the dictionary `L`, for instance `L["Vin07_L1"]` contains one matrix of Vinberg.



