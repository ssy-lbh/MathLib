cmake -S . -B build
cmake --build build
ctest --test-dir build -C Debug