cmake_minimum_required(VERSION 3.0)
project(NonHomog)




# File will default to direct sampling, continuous linear model, pure absorber solution

#########################################
######## USER SELECTED VARIABLES ########
#########################################
# Type of sampling method (default is direct)
option(REJECTION "Sample via thinning" OFF)

# Calculation type (default is linear)
option(QUADRATIC "Quadratic Model" ON)

# Function type (default is continuous)
option(PIECEWISE "Piecewise Constant Model" ON)

# Solution type (ONLY CHOOSE ONE) (geometry timing is prioritized - default is pure absorber solution)
option(S2 "Numerical Two-Ordinate SN Transport Method" OFF)
option(GEOMETRYTIME "Count CPU clock cycles required for certain number of geometry generations" OFF)
option(ANALYTICAL "One-dimensional analytical solution to transport equation - Requires BLAS/LAPACK" ON)
option(VOLFRAC "Set importance to one to compute volume fractions of generated materials" OFF)
#########################################
#########################################
#########################################




# Require standard math library
set(LIBRARIES_FOR_LINKING m)

# File will default to linear model with pure absorber solution
option(LINEAR "Linear Model" ON)
option(PUREABS "Numerical Pure Absorber Solution" ON)
# No need to change this ^^^

# Mutually exclusive variables
if(QUADRATIC)
    set(LINEAR OFF)
endif()

if(GEOMETRYTIME)
    set(PUREABS OFF)
    set(ANALYTICAL OFF)
    set(S2 OFF)
    set(VOLFRAC OFF)
endif()

if(S2)
    set(PUREABS OFF)
    set(ANALYTICAL OFF)
    set(VOLFRAC OFF)
endif()

if(ANALYTICAL)
    set(PUREABS OFF)
    set(VOLFRAC OFF)
    find_package(BLAS)
    find_package(LAPACK)
    if(LAPACK_FOUND AND BLAS_FOUND)
        set(LIBRARIES_FOR_LINKING ${LIBRARIES_FOR_LINKING} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
endif()

if(VOLFRAC)
    set(PUREABS OFF)
endif()

# Quadratic sampling is only via the thinning method
if(QUADRATIC AND NOT PIECEWISE)
    set(REJECTION OFF)
endif()

configure_file(
    "config.h.in"
    "config.h"
    @ONLY
)

include_directories(
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_SOURCE_DIR}"
)

add_compile_options(-O3 -Wall -Wextra -pedantic)

set(SRCS driver.c)

add_executable(NonHomog "${SRCS}")

target_link_libraries(NonHomog "${LIBRARIES_FOR_LINKING}")
