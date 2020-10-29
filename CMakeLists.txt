cmake_minimum_required( VERSION 3.0 )

project( cubature )

add_library( cubature SHARED 
    hcubature.c
    pcubature.c)
target_include_directories( cubature PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:.>)
  
add_executable( htest test.c )
target_link_libraries( htest cubature m )

add_executable( ptest test.c )
target_link_libraries( ptest cubature m )
target_compile_definitions( ptest PRIVATE PCUBATURE=1 )

include(GNUInstallDirs)
install( TARGETS cubature DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install( FILES cubature.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

