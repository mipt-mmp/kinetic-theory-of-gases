
add_library(phys STATIC 
universe.cpp
materialpoint.cpp
geometry.hpp materialpoint.cpp materialpoint.hpp physconstants.hpp units.hpp universe.cpp universe.hpp
real.hpp
)

target_include_directories(phys INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
