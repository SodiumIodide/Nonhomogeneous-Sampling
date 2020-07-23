cmake_minimum_required(VERSION 3.0)
project(NonHomog)




# File will default to direct sampling, continuous linear model, pure absorber solution

#########################################
######## USER SELECTED VARIABLES ########
#########################################
# Type of sampling method (default is direct)
option(REJECTION "Sample via thinning" OFF)

# Calculation type (default is linear)
option(QUADRATIC "Quadratic Model" OFF)

# Function type (default is continuous)
option(PIECEWISE "Piecewise Constant Model" OFF)

# Cox process (overrides previous selections)
option(COX_UNIFORM "Sampling via Doubly Stochastic (Cox) Process" OFF)
option(COX_EXPONENTIAL "Sampling via Cox Process from Exponential Distribution" OFF)
option(COX_GAUSSIAN "Sampling via Cox Process from Gaussian Distribution" ON)

# Constant chord length geometry generation (overrides previous selections)
option(CONSTANT "Constant chord length" OFF)

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
else()
    find_package(GSL)
    if (GSL_FOUND)
        set(LIBRARIES_FOR_LINKING ${LIBRARIES_FOR_LINKING} ${GSL_LIBRARIES})
    endif()
endif()

if(S2)
    set(PUREABS OFF)
    set(ANALYTICAL OFF)
    set(VOLFRAC OFF)
endif()

if(ANALYTICAL)
    set(PUREABS OFF)
    set(VOLFRAC OFF)
    find_package(LAPACK)
    if (LAPACK_FOUND)
        set(LIBRARIES_FOR_LINKING ${LIBRARIES_FOR_LINKING} ${LAPACK_LIBRARIES})
    endif()
endif()

if(VOLFRAC)
    set(PUREABS OFF)
endif()

# Quadratic sampling is only via the thinning method
if(QUADRATIC AND NOT PIECEWISE)
    set(REJECTION OFF)
endif()

if(COX_UNIFORM OR COX_EXPONENTIAL OR COX_GAUSSIAN)
    set(REJECTION OFF)
    set(QUADRATIC OFF)
    set(PIECEWISE OFF)
endif()

if (COX_UNIFORM)
    set(COX_EXPONENTIAL OFF)
    set(COX_GAUSSIAN OFF)
endif()

if (COX_EXPONENTIAL)
    set(COX_GAUSSIAN OFF)
endif()

if (CONSTANT)
    set(LINEAR OFF)
    set(QUADRATIC OFF)
    set(COX_UNIFORM OFF)
    set(COX_EXPONENTIAL OFF)
    set(COX_GAUSSIAN OFF)
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
