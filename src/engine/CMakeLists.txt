find_package(Qt6 REQUIRED COMPONENTS Concurrent)

add_library(phys STATIC 
chamber.cpp
geometry.hpp gasAtom.cpp gasAtom.hpp physconstants.hpp units.hpp chamber.cpp chamber.hpp ballsCollection.hpp ballsCollection.cpp
real.hpp
)

target_include_directories(phys INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(phys PRIVATE Qt6::Concurrent)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
