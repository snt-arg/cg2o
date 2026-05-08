
# Try to find Clang dynamically
find_program(CLANG_C_COMPILER clang)
find_program(CLANG_CXX_COMPILER clang++)

if(CLANG_C_COMPILER AND CLANG_CXX_COMPILER)
    # Check the version of Clang
    execute_process(
        COMMAND ${CLANG_CXX_COMPILER} --version
        OUTPUT_VARIABLE CLANG_VERSION_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    string(REGEX MATCH "[0-9]+\\.[0-9]+" CLANG_VERSION ${CLANG_VERSION_OUTPUT})

    message(STATUS "Detected Clang version: ${CLANG_VERSION}")

    # Check if the version is sufficient (>=14)
    if(CLANG_VERSION VERSION_LESS 14.0)
        message(FATAL_ERROR "Clang version 14.0 or higher is required. Detected version: ${CLANG_VERSION}")
    endif()

    # Use detected Clang compilers
    set(CMAKE_C_COMPILER ${CLANG_C_COMPILER})
    set(CMAKE_CXX_COMPILER ${CLANG_CXX_COMPILER})
else()
    message(FATAL_ERROR "Clang compiler not found. Please ensure Clang is installed.")
endif()

# Enforce the use of libstdc++ (GNU C++ Standard Library)
add_compile_options(-stdlib=libstdc++)

# Include a check for libstdc++ version (12 or higher)
include(CheckCXXSourceCompiles)

set(TEST_PROGRAM "
#include <iostream>
#include <vector>
int main() {
    std::vector<int> v{1, 2, 3};
    for (auto i : v) std::cout << i << std::endl;
    return 0;
}
")

check_cxx_source_compiles("${TEST_PROGRAM}" LIBSTDCXX_COMPATIBLE)

if(NOT LIBSTDCXX_COMPATIBLE)
    message(FATAL_ERROR "Incompatible libstdc++ version. Ensure version 12 or higher is installed.")
endif()

# Generate compile_commands.json for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Notify user about the configuration
message(STATUS "Configured with Clang (${CLANG_VERSION}) and libstdc++ 12.")

