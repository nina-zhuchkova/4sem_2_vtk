cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11")

project(vtkdemo)
add_executable(cubic2d cubic2d.cpp)
add_executable(cubic3d cubic3d.cpp)
add_executable(tetr3d tetr3d.cpp)
add_executable(tetr3d_from_msh tetr3d_from_msh.cpp)
add_executable(bug3d_wave bug3d_wave.cpp)
add_executable(bug3d_wave_move bug3d_wave_move.cpp)

# пути к библиотекам
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)
set(GMSH_LIBRARIES gmsh)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# ... тут только VTK
target_link_libraries(cubic2d ${VTK_LIBRARIES})
target_link_libraries(cubic3d ${VTK_LIBRARIES})

# ... а тут ещё и gmsh
target_link_libraries(tetr3d ${VTK_LIBRARIES})
target_link_libraries(tetr3d ${GMSH_LIBRARIES})
target_link_libraries(tetr3d_from_msh ${VTK_LIBRARIES})
target_link_libraries(tetr3d_from_msh ${GMSH_LIBRARIES})
target_link_libraries(bug3d_wave ${VTK_LIBRARIES})
target_link_libraries(bug3d_wave ${GMSH_LIBRARIES})
target_link_libraries(bug3d_wave_move ${VTK_LIBRARIES})
target_link_libraries(bug3d_wave_move ${GMSH_LIBRARIES})
