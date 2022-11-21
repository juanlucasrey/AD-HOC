#!/bin/zsh
cd build
./tests/adhoc_test
/usr/local/Cellar/llvm/15.0.0/bin/llvm-profdata merge -sparse default.profraw -o default.profdata
/usr/local/Cellar/llvm/15.0.0/bin/llvm-cov show ./tests/adhoc_test -instr-profile=default.profdata -format=html >> cov.html

/usr/local/Cellar/llvm/15.0.0/bin/llvm-cov report ./tests/adhoc_test -instr-profile=default.profdata -show-instantiation-summary