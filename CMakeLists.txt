cmake_minimum_required(VERSION 3.0)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
# Not strictly necessary, but will create executables that are
# aware of ALPSCore location
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

project(alps_cthyb)

# Disable in-source builds
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Request the ALPSCore package:
find_package(ALPSCore COMPONENTS accumulators params mc hdf5 gf)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI)
find_package(NFFT3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS -DNDEBUG")

file(GLOB_RECURSE HEADER_DIRS "src/*")  # src以下のすべてのディレクトリを収集
foreach(HEADER_DIR ${HEADER_DIRS})
    if(IS_DIRECTORY ${HEADER_DIR})
        include_directories(${HEADER_DIR})  # ディレクトリをインクルードパスに追加
    endif()
endforeach()

include_directories(
    ${HEADER_FILES}
    ${CMAKE_BINARY_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${NFFT3_INCLUDE_DIRS}
)

file(GLOB_RECURSE ALPS_CTHYB_SOURCE "src/**/*.cpp")

add_library(alps_hyb_exp ${ALPS_CTHYB_SOURCE})
target_link_libraries(alps_hyb_exp ${ALPSCore_LIBRARIES} ${MPI_CXX_LIBRARIES} ${NFFT3_LIBRARIES})
set_target_properties(alps_hyb_exp PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(EXEC_NAME alps_cthyb)
add_executable(${EXEC_NAME} src/hybmain.cpp)
target_link_libraries(${EXEC_NAME} alps_hyb_exp ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES})
set_target_properties(${EXEC_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
install (TARGETS ${EXEC_NAME} DESTINATION bin )

### Additional installation files
# configuration file for module managers, such as LMOD
configure_file("${CMAKE_SOURCE_DIR}/misc/alps_cthyb.lmod.in" "${CMAKE_BINARY_DIR}/alps_cthyb.lmod")
