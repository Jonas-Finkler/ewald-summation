# Copyright (C) 2020 Jonas A. Finkler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


cmake_minimum_required(VERSION 3.6)
project(ewald-summation Fortran)
enable_language(Fortran)


option(INTEL "compile with ifort and mkl" OFF)
option(OPENMP "compile with openmp support" OFF)
option(DEBUG "compile with debug flags" OFF)

if(INTEL)
    set(CMAKE_Fortran_COMPILER ifort)
else()
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

if(INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin")
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -g -traceback -O0 -fpe0") # -Ofast -ipo
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -mkl -Ofast") # -ipo
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")# -parallel")
    endif()
else()
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall -Wconversion -g -fbacktrace -ffpe-trap=zero,invalid")
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -march=native -flto -fno-strict-aliasing")
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()



add_executable(
        ewald-summation.x

        src/testEwaldSummation.f90
        src/constants/constants.f90
        src/util/util.f90
        src/configuration/parameters.f90
        src/ewaldSummation/ewaldSummation.f90
        src/atoms/pointsets.f90
        src/precision.f90
        src/math/linalg.f90
        src/math/random.f90
)

target_link_libraries(ewald-summation.x ${libraries})
