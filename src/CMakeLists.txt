cmake_minimum_required( VERSION 2.8.2 )

# include_directories( ${FLEXBAR_SOURCE_DIR}/include )
# link_directories( ${FLEXBAR_SOURCE_DIR}/lib )

if( CMAKE_SIZEOF_VOID_P MATCHES "8" )
	message( STATUS "Flexbar 64 bit architecture" )
else()
	message( STATUS "Flexbar 32 bit architecture" )
endif()

add_executable( flexbar Flexbar.cpp )
target_link_libraries( flexbar tbb )

find_package( ZLIB )
if( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( flexbar ${ZLIB_LIBRARIES} )
	add_definitions( "-DSEQAN_HAS_ZLIB=1" )
else()
	message( STATUS "Build will not support zlib!" )
endif()

find_package( BZip2 )
if( BZIP2_FOUND )
	include_directories( ${BZIP2_INCLUDE_DIRS} )
	target_link_libraries( flexbar ${BZIP2_LIBRARIES} )
	add_definitions( "-DSEQAN_HAS_BZIP2=1" )
else()
	message( STATUS "Build will not support bzip2!" )
endif()

