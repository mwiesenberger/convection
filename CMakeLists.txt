cmake_minimum_required(VERSION 2.8.8)
project(Convection)
##################set up Qt cmake environment#####################
set(convection_SOURCES main.cpp mainwindow.cpp window.cpp)
set(convection_HEADERS mainwindow.h window.h) # all headers for moc compiler
set(convection_FORMS mainwindow.ui window.ui) # all headers for moc compiler
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set( CMAKE_AUTOMOC ON)
#set( CMAKE_AUTOUIC ON)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
#QT5_WRAP_CPP(convection_HEADERS_MOC ${convection_HEADERS})
QT5_WRAP_UI(convection_FORMS_HEADERS ${convection_FORMS})
#include(${QT_USE_FILE})
#add_definitions(${QT_DEFINITIONS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
#
add_executable(convection
 ${convection_SOURCES}
 ${convection_FORMS_HEADERS})

qt5_use_modules( convection Widgets Core Gui)
########################set up the rest#######################
set( CMAKE_BUILD_TYPE Debug)
find_package( OpenCV REQUIRED)  # find well-known library for which cmake has config files for many platforms
if( WIN32)
find_library( FFTW_LIB NAMES fftw3-3 fftw3 PATHS C:/MinGW/msys/1.0/local/bin NO_DEFAULT_PATH)   # find system library
elseif( UNIX)
find_library( FFTW_LIB NAMES fftw3)
endif()
find_package( OpenMP )

if( UNIX)
include_directories( /home/matthias/include)
include_directories( /home/matthias/Projekte/feltor/inc/)
elseif( WIN32)
INCLUDE_DIRECTORIES( C:/Users/Matthias/Documents/GitHub)
INCLUDE_DIRECTORIES( C:/Users/Matthias/Documents/GitHub/feltor/inc)
INCLUDE_DIRECTORIES( C:/MinGW/msys/1.0/local/include)
endif()

target_link_libraries(convection ${OpenCV_LIBS})
target_link_libraries(convection ${QT_LIBRARIES})
target_link_libraries(convection ${FFTW_LIB})

set( CMAKE_CXX_FLAGS "-std=c++0x")
###############Set compiler flags for OpenMP support##############
#if( OpenMP)
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()
