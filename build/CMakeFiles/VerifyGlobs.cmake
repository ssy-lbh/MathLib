# CMAKE generated file: DO NOT EDIT!
# Generated by CMake Version 3.22
cmake_policy(SET CMP0009 NEW)

# math_lib_src at lib/CMakeLists.txt:3 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "D:/Code/C++/MathLib/lib/*.cpp")
set(OLD_GLOB
  "D:/Code/C++/MathLib/lib/math/cipolla.cpp"
  "D:/Code/C++/MathLib/lib/math/crt.cpp"
  "D:/Code/C++/MathLib/lib/math/index_calculus.cpp"
  "D:/Code/C++/MathLib/lib/math/math_base.cpp"
  "D:/Code/C++/MathLib/lib/math/miller_rabin.cpp"
  "D:/Code/C++/MathLib/lib/math/pohlig_hellman.cpp"
  "D:/Code/C++/MathLib/lib/math/pollard_rho.cpp"
  "D:/Code/C++/MathLib/lib/math/prime_count.cpp"
  "D:/Code/C++/MathLib/lib/math/primitive_root.cpp"
  "D:/Code/C++/MathLib/lib/math/sieves.cpp"
  "D:/Code/C++/MathLib/lib/math/totient.cpp"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "D:/Code/C++/MathLib/build/CMakeFiles/cmake.verify_globs")
endif()

# math_lib_src at lib/CMakeLists.txt:3 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "D:/Code/C++/MathLib/lib/*c")
set(OLD_GLOB
  "D:/Code/C++/MathLib/lib/hash/sha1.c"
  "D:/Code/C++/MathLib/lib/hash/sha256.c"
  "D:/Code/C++/MathLib/lib/string/kmp.c"
  "D:/Code/C++/MathLib/lib/util/zipmap.c"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "D:/Code/C++/MathLib/build/CMakeFiles/cmake.verify_globs")
endif()
