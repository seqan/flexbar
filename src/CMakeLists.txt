cmake_minimum_required( VERSION 2.8.2 )

if( CMAKE_SIZEOF_VOID_P MATCHES "8" )
	message( STATUS "Flexbar 64 bit architecture" )
	# if( WIN32 )
	# 	link_directories( ${FLEXBAR_SOURCE_DIR}/lib/win64 )
	# elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	# 	link_directories( ${FLEXBAR_SOURCE_DIR}/lib/macosx )
	# else() # assuming linux, or adjust to your system's tbb lib
	# 	link_directories( ${FLEXBAR_SOURCE_DIR}/lib/linux64 )
	# endif()
else()
	message( STATUS "Flexbar 32 bit architecture" )
	# if( WIN32 )
	# 	link_directories( ${FLEXBAR_SOURCE_DIR}/lib/win32 )
	# else()
	# 	message( FATAL_ERROR "Platform not supported!" )
	# endif()
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


set( SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} flexbar CACHE INTERNAL "" )

# set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/Users/jtr/Downloads/seqan-trunk/util/cmake" )
# set( SEQAN_FIND_DEPENDENCIES ZLIB BZip2 )
# find_package( SeqAn REQUIRED )

# add_definitions( ${SEQAN_DEFINITIONS} )
# include_directories( ${SEQAN_INCLUDE_DIRS} )
# target_link_libraries( flexbar ${SEQAN_LIBRARIES} )
# set( SEQAN_HAS_ZLIB FALSE )

# find_package( TBB REQUIRED )
# if( NOT TBB_FOUND )
# 	message( FATAL_ERROR "TBB library not found." )
# endif()
