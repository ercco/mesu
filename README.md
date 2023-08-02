# Subnetwork enumeration algorithms for multilayer networks

See https://doi.org/10.48550/arXiv.2308.00083 for the theoretical background.

This repository contains two algorithms for enumerating subnetworks of multilayer networks, both implemented in Python and in C++. For the Python implementation, the [pymnet](http://www.mkivela.com/pymnet/) library is used. C++ version uses Boost Graph Library.

Compile c++ with:
`g++ -DN_ASPECTS=n_aspects -std=c++17 -O3 mesu.cpp -o mesu_n_aspects.out`
where **n_aspects** should be replaced with the desired number of aspects.

Or:
`make`
will compile versions from 1 to 6 aspects.

Use c++ with:
`mesu_d.out inputfile outputfile 'subnet_size_in_aspect_0,subnet_size_in_aspect_1,...,subnet_size_in_aspect_d'`

To run c++ benchmarks with `run_benchmark_models_cpp` (in `benchmarks.py`) there need to be compiled versions of `mesu.cpp` with names mesu\_**n_aspects**.out.
