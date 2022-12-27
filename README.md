# MESU
Multilayer EnumerateSubgraphs

Compile c++ with:
`g++ -DN_ASPECTS=n_aspects -std=c++17 -O3 mesu.cpp -o mesu_n_aspects.out`
where n_aspects should be replaced with the desired number of aspects.

Or:
`make`
will compile versions from 1 to 6 aspects.

Use c++ with:
`mesu.out inputfile outputfile 'subnet_size_in_aspect_0,subnet_size_in_aspect_1,...,subnet_size_in_aspect_d'`

To run c++ benchmarks with `run_benchmark_models_cpp` (in `benchmarks.py`) there need to be compiled versions of `mesu.cpp` with names **mesu_N.out** where **N** corresponds to `#define N_ASPECTS N` constant in `mesu.cpp`.
