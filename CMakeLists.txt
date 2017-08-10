cmake_minimum_required( VERSION 3.0 )

project( cubature )

add_library( cubature SHARED 
    hcubature.c
    pcubature.c)

add_executable( htest test.c )
target_link_libraries( htest cubature m )

add_executable( ptest test.c )
target_link_libraries( ptest cubature m )
target_compile_definitions( ptest PRIVATE PCUBATURE=1 )

include(GNUInstallDirs)
install( TARGETS cubature DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install( FILES cubature.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

