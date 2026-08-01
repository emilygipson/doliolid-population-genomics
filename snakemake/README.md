# Nuclear lcWGS workflow

ANGSD steps for the *D. gegenbauri* nuclear low-coverage data, n=118. SAF, folded
SFS, thetas, pairwise Fst, PCAngsd. Mapping and dedup happen upstream in
`scripts/bash/` and aren't part of this.

Full pipeline writeup: `markdowns/nuclear_lcwgs_pipeline.md`.

## Groups and pairs

Groups are bamlists, named in `config.yaml`. Every per-group rule runs across all
of them. Pairs name two groups and get Fst under both `-whichFst 0` (Reynolds)
and `-whichFst 1` (Bhatia).
