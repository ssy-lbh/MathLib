project(math_lib)

file(GLOB_RECURSE math_lib_src CONFIGURE_DEPENDS *c *.cpp)

add_library(${PROJECT_NAME} STATIC ${math_lib_src})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
