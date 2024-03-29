# - Try to find MUMPS
# Once done this will define
#  MUMPS_FOUND - System has Mumps
#  MUMPS_INCLUDE_DIRS - The Mumps include directories
#  MUMPS_LIBRARY_DIRS - The library directories needed to use Mumps
#  MUMPS_LIBRARIES    - The libraries needed to use Mumps

if (MUMPS_INCLUDE_DIR)
  # in cache already
  SET(MUMPS_FIND_QUIETLY TRUE)
endif (MUMPS_INCLUDE_DIR)

find_path(MUMPS_INCLUDE_DIR
        NAMES dmumps_c.h
        HINTS "/home/detu/coding/ThirdParty-Mumps/MUMPS/include"
        "/usr/local/include/coin-or/mumps"
        )

find_library( MUMPS_LIBRARY_DEBUG
        coinmumps
        HINTS "/home/detu/coding/ThirdParty-Mumps/MUMPS/.libs"
        "/usr/lib"
        )


include(SelectLibraryConfigurations)
select_library_configurations( MUMPS )

set(MUMPS_INCLUDE_DIRS "${MUMPS_INCLUDE_DIR}" )
set(MUMPS_LIBRARIES "${MUMPS_LIBRARY}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MUMPS  DEFAULT_MSG
        MUMPS_LIBRARY MUMPS_INCLUDE_DIR)

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY )