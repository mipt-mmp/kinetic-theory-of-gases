find_package(Qt6 REQUIRED COMPONENTS Core Widgets)
qt_standard_project_setup()

qt_add_library(phys STATIC
    chamber.cpp
    geometry.hpp
    gasAtom.cpp
    gasAtom.hpp
    physconstants.hpp
    units.hpp
    chamber.cpp
    chamber.hpp
    real.hpp
)

target_link_libraries(phys PRIVATE Qt6::Core Qt6::Widgets)

target_include_directories(phys INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
