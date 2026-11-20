cmake_minimum_required(VERSION 3.20)
project(day001h_asian_options LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

find_package(Threads REQUIRED)
find_library(MATH_LIB m)

set(ASIAN_C_FLAGS -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)

add_library(asian_c_lib
  src/options.c
  src/svg.c
)
target_include_directories(asian_c_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(asian_c_lib PRIVATE ${ASIAN_C_FLAGS})
target_link_libraries(asian_c_lib PUBLIC Threads::Threads)
if(MATH_LIB)
  target_link_libraries(asian_c_lib PUBLIC ${MATH_LIB})
endif()

if(ENABLE_UBSAN)
  target_compile_options(asian_c_lib PUBLIC -fsanitize=undefined -fno-omit-frame-pointer)
  target_link_options(asian_c_lib PUBLIC -fsanitize=undefined)
endif()

add_executable(day001h_run app/main.c)
target_link_libraries(day001h_run PRIVATE asian_c_lib)
target_compile_options(day001h_run PRIVATE ${ASIAN_C_FLAGS})

enable_testing()
add_executable(day001h_tests tests/test_options.c)
target_link_libraries(day001h_tests PRIVATE asian_c_lib)
target_compile_options(day001h_tests PRIVATE ${ASIAN_C_FLAGS})
add_test(NAME asian_c_tests COMMAND day001h_tests)
