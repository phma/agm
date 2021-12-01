project(agm)
cmake_minimum_required(VERSION 3.8.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

add_executable(
  agm main.cpp agm.cpp angle.cpp ps.cpp ldecimal.cpp
)

# Define NO_INSTALL when compiling for fuzzing. This avoids the error
# "The install of the agm target requires changing an RPATH", which
# occurs when using the AFL compiler wrapper with the Ninja generator.
# There is no need to install a binary built for fuzzing.
if (NOT DEFINED NO_INSTALL)
install(TARGETS agm DESTINATION bin)
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(AGM_MAJOR_VERSION 0)
set(AGM_MINOR_VERSION 1)
set(AGM_PATCH_VERSION 0)
set(AGM_VERSION ${AGM_MAJOR_VERSION}.${AGM_MINOR_VERSION}.${AGM_PATCH_VERSION})
set(COPY_YEAR 2021)

include_directories(${PROJECT_BINARY_DIR})
configure_file(config.h.in config.h)

set(CPACK_PACKAGE_VERSION_MAJOR ${AGM_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${AGM_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${AGM_PATCH_VERSION})
set(CPACK_SOURCE_IGNORE_FILES /\\\\.git;.*~)
include(CPack)
