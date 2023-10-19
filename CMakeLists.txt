cmake_minimum_required(VERSION 3.20)

project(celestial-mechanics VERSION 0.1.0 LANGUAGES CXX)

include(cmake/Logging.cmake)
include(cmake/Ccache.cmake)
include(cmake/Helpers.cmake)
include(cmake/CheckCompiler.cmake)

add_subdirectory(lib)

include(cmake/CompileOptions.cmake)
#include(cmake/Sanitize.cmake)
include(cmake/Development.cmake)

# Must be after lib not to affect them.
# include(cmake/Werror.cmake)

add_subdirectory(src)
