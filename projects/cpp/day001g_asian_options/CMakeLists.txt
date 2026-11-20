cmake_minimum_required(VERSION 3.20)
project(day001g_asian_options LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

add_library(asian_lib
  src/options.cpp
  src/svg.cpp
)
target_include_directories(asian_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(asian_lib PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(ENABLE_SANITIZERS)
  target_compile_options(asian_lib PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(asian_lib PUBLIC -fsanitize=address,undefined)
endif()

add_executable(day001g_run app/main.cpp)
target_link_libraries(day001g_run PRIVATE asian_lib)
target_compile_options(day001g_run PRIVATE -Wall -Wextra -Wpedantic -Werror)

include(CTest)
include(FetchContent)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Policy version for FetchContent deps")
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG v2.4.11
)
FetchContent_MakeAvailable(doctest)

add_executable(day001g_tests tests/test_options.cpp)
target_link_libraries(day001g_tests PRIVATE asian_lib doctest::doctest)
target_compile_options(day001g_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_test(NAME asian_tests COMMAND day001g_tests)
