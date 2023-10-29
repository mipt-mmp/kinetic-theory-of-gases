find_package(Qt6 REQUIRED COMPONENTS Core Widgets)
qt_standard_project_setup()

qt_add_library(phys STATIC
    gasAtom.cpp
    gasAtom.hpp
    chamber.cpp
    chamber.hpp
    sync.cpp
    sync.hpp
    units.hpp
    geometry.hpp
    physconstants.hpp
    real.hpp
)

target_link_libraries(phys PRIVATE Qt6::Core Qt6::Widgets)

target_include_directories(phys INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
