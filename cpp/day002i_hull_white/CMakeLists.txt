cmake_minimum_required(VERSION 3.20)
project(day002i_hull_white LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

add_library(hw_lib
  src/hull_white.cpp
  src/trinomial_tree.cpp
  src/calibration.cpp
  src/bermudan.cpp
  src/montecarlo.cpp
  src/svg.cpp
)
target_include_directories(hw_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(hw_lib PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(ENABLE_SANITIZERS)
  target_compile_options(hw_lib PUBLIC -fsanitize=undefined -fno-omit-frame-pointer)
  target_link_options(hw_lib PUBLIC -fsanitize=undefined)
endif()

add_executable(day002i_run app/main.cpp)
target_link_libraries(day002i_run PRIVATE hw_lib)
target_compile_options(day002i_run PRIVATE -Wall -Wextra -Wpedantic -Werror)

include(CTest)
include(FetchContent)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Policy version for FetchContent deps")
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG v2.4.11
)
FetchContent_MakeAvailable(doctest)

add_executable(day002i_tests tests/test_hw.cpp)
target_link_libraries(day002i_tests PRIVATE hw_lib doctest::doctest)
target_compile_options(day002i_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_test(NAME hw_tests COMMAND day002i_tests)
