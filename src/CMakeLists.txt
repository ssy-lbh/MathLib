project(math_lib)

file(GLOB_RECURSE math_lib_src CONFIGURE_DEPENDS *c *.cpp)

add_library(${PROJECT_NAME} STATIC ${math_lib_src})

target_link_directories(${PROJECT_NAME} PUBLIC ${LIBRARY_PATH})

target_link_libraries(${PROJECT_NAME} PUBLIC gmp::gmp gmp-10 mpfr-4)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
