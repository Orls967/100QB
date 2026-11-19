cmake_minimum_required(VERSION 3.20)
project(day002h_hawkes_orderbook LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

add_library(hawkes_ob_lib
  src/hawkes.cpp
  src/orderbook.cpp
  src/market_sim.cpp
  src/svg.cpp
)
target_include_directories(hawkes_ob_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(hawkes_ob_lib PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(ENABLE_SANITIZERS)
  target_compile_options(hawkes_ob_lib PUBLIC -fsanitize=undefined -fno-omit-frame-pointer)
  target_link_options(hawkes_ob_lib PUBLIC -fsanitize=undefined)
endif()

add_executable(day002h_run app/main.cpp)
target_link_libraries(day002h_run PRIVATE hawkes_ob_lib)
target_compile_options(day002h_run PRIVATE -Wall -Wextra -Wpedantic -Werror)

include(CTest)
include(FetchContent)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Policy version for FetchContent deps")
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG v2.4.11
)
FetchContent_MakeAvailable(doctest)

add_executable(day002h_tests tests/test_hawkes_ob.cpp)
target_link_libraries(day002h_tests PRIVATE hawkes_ob_lib doctest::doctest)
target_compile_options(day002h_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_test(NAME hawkes_ob_tests COMMAND day002h_tests)
