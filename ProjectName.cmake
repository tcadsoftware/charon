SET(PROJECT_NAME tcad-charon)
SET(${PROJECT_NAME}_GENERATE_REPO_VERSION_FILE TRUE CACHE BOOL
  "Generate repo version file to tie repos together at certain commits.  "
  "Set in ProjectName.cmake.")
SET(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT ON)
SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE Nightly CACHE STRING
  "Enable the extra repositories.  Set in ProjectName.cmake.")
SET(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES TRUE CACHE BOOL
  "Ignore extra repositories that are not present.  Set in ProjectName.cmake.")
SET(${PROJECT_NAME}_ENABLE_GCOVR FALSE CACHE BOOL
  "Turn on GCOVR support.")
