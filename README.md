# MESU
Multilayer EnumerateSubgraphs

Compile c++ with:
`g++ -std=c++17 -O3 mesu.cpp -o mesu.out`

Use c++ with:
`mesu.out inputfile outputfile 'subnet_size_in_aspect_1,subnet_size_in_aspect_2,...,subnet_size_in_aspect_d'`

To run c++ benchmarks with `run_benchmark_models_cpp` (in `benchmarks.py`) there need to be compiled versions of `mesu.cpp` with names **mesu_N.out** where **N** corresponds to `#define N_ASPECTS N` constant in `mesu.cpp`.
