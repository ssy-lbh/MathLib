project(math_test)

function(add_math_test name)
	add_executable(${name} ${name}.cpp)
	target_include_directories(${name} PUBLIC ${SOURCE_PATH})
	target_link_libraries(${name} PUBLIC math_lib)
	add_test(NAME ${name} COMMAND ${name})
endfunction(add_math_test)

add_math_test(mod)
add_math_test(algo)
add_math_test(prime)
add_math_test(bignum)
add_math_test(poly)
