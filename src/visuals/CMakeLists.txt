
find_package(Qt6 REQUIRED COMPONENTS Core Widgets)
qt_standard_project_setup()

qt_add_executable(mkt
    main.cpp
    mainwindow.cpp
    mainwindow.hpp
    mainwindow.cpp mainwindow.ui
    chamberdisplayer.hpp chamberdisplayer.cpp
    physicsthread.hpp physicsthread.cpp
)

target_link_libraries(mkt PRIVATE phys Qt6::Core Qt6::Widgets)
