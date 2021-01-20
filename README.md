# 𝐕𝐢𝐧𝐧𝐲 𝐭𝐡𝐞 𝐜𝐨𝐱𝐛𝐨𝐲 & 𝐉𝐮𝐥𝐢𝐚

### Learning to herd Coxeter diagrams in the Hyperbolic Plains since 2021

***

# Port of [aperep/vinberg-algorithm](https://github.com/aperep/vinberg-algorithm) to julia, plus some changes


## Still being tested: probably many bugs left!


* Finite volume check for Coxeter diagrams implemented “from scratch” using Guglielmetti's thesis as theoretical resource and [`CoxIter`](https://github.com/rgugliel/CoxIter) as comparison point (plus a majority of tests from there)
* The main bulk of Vinberg's algorithm (and associated procedures) ported from [B&P](https://github.com/aperep/vinberg-algorithm)
* Using https://raw.githubusercontent.com/chethega/StaticArrays.jl/fb0350012f01db4021d60906357e949333ec5d93/src/SBitSet.jl for a fixed size bitset implementation. **TODO** contact author to ask if we can use their code, and under what licence (see [here](https://github.com/chethega/StaticBitsets.jl/issues/1)).
