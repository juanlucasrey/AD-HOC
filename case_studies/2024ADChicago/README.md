# AD 2024 Conference Paper

These instructions will build small programs differentiating the Black-Scholes formula at different orders, using different packages. Assuming you want to run 1.000.000 different Black-Scholes prices and calculate its derivatives with a program "package_orderX", make sure you run from the terminal with the following command:\
ITERATIONS=1000000 ./package_orderX

For example if you want to run with 1.000.000 Black-Scholes with CoDiPack at order 2, run the following command:\
ITERATIONS=1000000 ./codi_order2

You will get as a result the average price and derivatives, with the execution time in ms as follows:\
7.7723986486630947\
...\
-1.4731780924569864\
iterations: 1000000\
adhoc order 2 time (ms): 500

assuming /path/g++ is your path to GCC++, build the timing executables with the following commands:

## Primal calculation (no derivatives)
compile with:

/path/g++ -o primal primal.cpp -std=c++20 -O3

## AD-HOC
compile with:

/path/g++ -o adhoc_order1 adhoc_order1.cpp -std=c++20 -O3 -I../../include\
/path/g++ -o adhoc_order2 adhoc_order2.cpp -std=c++20 -O3 -I../../include\
/path/g++ -o adhoc_order3 adhoc_order3.cpp -std=c++20 -O3 -I../../include\
/path/g++ -o adhoc_order4 adhoc_order4.cpp -std=c++20 -O3 -I../../include\
/path/g++ -o adhoc_order5 adhoc_order5.cpp -std=c++20 -O3 -I../../include

note: make sure you use GCC as clang will break for order 3 onwards. Also for order 5, at the time of writing, compilation will take a while. Current work is being done to reduce compilations times and make sure clang works for higher orders.

## [ADOL-C](https://github.com/coin-or/ADOL-C)
Assuming [ADOL-C](https://github.com/coin-or/ADOL-C) is installed in the folder /install/adolc, compile with:

/path/g++ -o adolc_order1 adolc_order1.cpp -std=c++20 -O3 -I/install/adolc/include -L/install/adolc/lib  -rpath /install/adolc/lib -ladolc\
/path/g++ -o adolc_order2 adolc_order2.cpp -std=c++20 -O3 -I/install/adolc/include -L/install/adolc/lib  -rpath /install/adolc/lib -ladolc\
/path/g++ -o adolc_order3 adolc_order3.cpp -std=c++20 -O3 -I/install/adolc/include -L/install/adolc/lib  -rpath /install/adolc/lib -ladolc\
/path/g++ -o adolc_order4 adolc_order4.cpp -std=c++20 -O3 -I/install/adolc/include -L/install/adolc/lib  -rpath /install/adolc/lib -ladolc\
/path/g++ -o adolc_order5 adolc_order5.cpp -std=c++20 -O3 -I/install/adolc/include -L/install/adolc/lib  -rpath /install/adolc/lib -ladolc

## [dco/c++](https://nag.com/automatic-differentiation)
Assuming [dco/c++](https://nag.com/automatic-differentiation/) is installed in the folder /install/dco, compile with:

/path/g++ -o dco_order1 dco_order1.cpp -std=c++20 -O3 -I/install/dco/include\
/path/g++ -o dco_order2 dco_order2.cpp -std=c++20 -O3 -I/install/dco/include

note: dco/c++ is a commercial package so you will need to purchase a licence.

## [dco/map](https://nag.com/automatic-differentiation)
Assuming [dco/map](https://nag.com/automatic-differentiation/) is installed in the folder /install/dcomap, compile with:

/path/g++ -o dcomap_order1 dco_map_order1.cpp -std=c++20 -O3 -I/install/dcomap/src

note: dco/map is a commercial package so you will need to purchase a licence.

## [CoDiPack](https://github.com/SciCompKL/CoDiPack)
Assuming [CoDiPack](https://github.com/SciCompKL/CoDiPack) is installed in the folder /install/codi, compile with:

/path/g++ -o codi_order1 codi_order1.cpp -std=c++20 -O3 -I/install/codi\
/path/g++ -o codi_order2 codi_order2.cpp -std=c++20 -O3 -I/install/codi

## [Enzyme](https://github.com/EnzymeAD/Enzyme)
Assuming [Enzyme](https://github.com/EnzymeAD/Enzyme) is installed in the folder /install/enzyme, compile with:

/path/clang++ -o enzyme_order1 enzyme_order1.cpp -O3 -Xclang -load -Xclang /install/enzyme/lib/ClangEnzyme-19.dylib -std=c++20

note: Enzyme requires clang++ instead of GCC

## [Tapenade](https://gitlab.inria.fr/tapenade/tapenade)
[Tapenade](https://gitlab.inria.fr/tapenade/tapenade) does not support c++ at the time of writing, so a different header, written in C ("black_scholes_c.h") has been used to generate the Tapenade code. These generated files are stored in the folder Tapenade. Compile with:

/path/g++ -o tapenade_order1 tapenade_order1.c tapenade/c_version_bs_b.c tapenade/adStack.c -O3\
/path/g++ -o tapenade_order2 tapenade_order2.cpp tapenade/c_version_bs_d_b_corrected.c tapenade/adStack.c -O3

## Python correctness verification
The correctness of the derivative calculations for the first iteration can be verified using a python script "sympy_value_check.py". Make sure you have python3 with sympy installed and run:\
python3 sympy_value_check.py

The values being printed should correspond to the values of any program being run, with ITERATIONS=1 (or the ITERATIONS argument being absent since 1 is the default value)
