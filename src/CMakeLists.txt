cmake_minimum_required( VERSION 2.8.0 )

# include_directories( ${FLEXBAR_SOURCE_DIR}/include )
# link_directories( ${FLEXBAR_SOURCE_DIR}/lib )

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


