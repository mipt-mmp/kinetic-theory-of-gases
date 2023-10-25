
add_library(phys STATIC 
chamber.cpp
geometry.hpp gasAtom.cpp gasAtom.hpp physconstants.hpp units.hpp chamber.cpp chamber.hpp
real.hpp
)

target_include_directories(phys INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
