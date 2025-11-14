
cmake_minimum_required(VERSION 3.16)

project(LinkTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(../XFoil-lib/)
include_directories(../fl5-lib/)
include_directories(../fl5-lib/api/)
include_directories(/usr/include/qt6/QtCore)
include_directories(/usr/include/qt6)


link_directories(/usr/local/lib/XFoil/)
link_directories(/usr/local/lib/fl5-lib/)
link_directories(/usr/lib64/)


if (WIN32)
else ()
endif (WIN32)

add_executable(LinkTest fl5link.cpp)

target_link_libraries(LinkTest XFoil fl5-lib Qt6Core)

include(GNUInstallDirs)
install(TARGETS LinkTest
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
