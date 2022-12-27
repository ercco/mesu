all: mesu.cpp
	for N_ASPECTS in 1 2 3 4 5 6 ; do g++ -DN_ASPECTS=$${N_ASPECTS} -std=c++17 -O3 mesu.cpp -o mesu_$${N_ASPECTS}.out ; done
