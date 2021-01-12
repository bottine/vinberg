# Port of [aperep/vinberg-algorithm](https://github.com/aperep/vinberg-algorithm) to julia, plus some changes


* Finite volume check for Coxeter diagrams implemented “from scratch” using Guglielmetti's thesis as theoretical resource and [`CoxIter`](https://github.com/rgugliel/CoxIter) as comparison point (plus a majority of tests from there)
* The main bulk of Vinberg's algorithm (and associated procedures) ported from [B&P](https://github.com/aperep/vinberg-algorithm)

