cmake_minimum_required(VERSION 2.8)
project(Convection)
##################set up Qt cmake environment#####################
find_package( Qt4 REQUIRED)

set(convection_SOURCES main.cpp mainwindow.cpp)
set(convection_HEADERS mainwindow.h) # all headers for moc compiler
set(convection_FORMS mainwindow.ui) # all headers for moc compiler

QT4_WRAP_CPP(convection_HEADERS_MOC ${convection_HEADERS})
QT4_WRAP_UI(convection_FORMS_HEADERS ${convection_FORMS})
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(convection
 ${convection_SOURCES}
 ${convection_HEADERS_MOC}
 ${convection_FORMS_HEADERS})

########################set up the rest#######################
set( CMAKE_BUILD_TYPE Release)
find_package( OpenCV REQUIRED)  # find well-known library for which cmake has config files for many platforms
find_library( FFTW3 fftw3)   # find system library
find_package( OpenMP REQUIRED )

include_directories( /home/matthias/include)
include_directories( /home/matthias/Projekte/feltor/inc/)

target_link_libraries(convection ${OpenCV_LIBS})
target_link_libraries(convection ${QT_LIBRARIES})
target_link_libraries(convection ${FFTW3})

set( CMAKE_CXX_FLAGS "-std=c++0x")
###############Set compiler flags for OpenMP support##############
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
