TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( BoostCharonLib
  REQUIRED_HEADERS boost/version.hpp boost/mpl/at.hpp
  MUST_FIND_ALL_HEADERS
  REQUIRED_LIBS_NAMES boost_filesystem boost_regex boost_mpi boost_exception boost_serialization
  MUST_FIND_ALL_LIBS
  )
