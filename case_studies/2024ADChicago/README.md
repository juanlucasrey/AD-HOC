## AD 2024 Conference Paper

assuming path/g++ is your path to GCC, build with the executables with the following commands:

for the primal calculation timing:
path/g++ -o primal primal.cpp -std=c++20 -O3

for AD-HOC timing:
path/g++ -o adhoc_order1 adhoc_order1.cpp -std=c++20 -O3
path/g++ -o adhoc_order2 adhoc_order2.cpp -std=c++20 -O3
path/g++ -o adhoc_order3 adhoc_order3.cpp -std=c++20 -O3
path/g++ -o adhoc_order4 adhoc_order4.cpp -std=c++20 -O3
path/g++ -o adhoc_order5 adhoc_order5.cpp -std=c++20 -O3
make sure you GCC as clang will break for order 3 onwards. Also for order 5, at the time of writing, compilation will take a while. Current work is being done to reduce compilations times.

for NAG-DCO timing, make sure NAG-DCO is installed in the folder ./packages/dco:
path/g++ -o nag_order1 nag_order1.cpp -std=c++20 -O3


for ADOL-C timing, make sure ADOL-C is installed in
note: ADOL-C does not support erfc at the time of writing so a special version without erfc has been written in the file "black_scholes_no_erfc.hpp"

