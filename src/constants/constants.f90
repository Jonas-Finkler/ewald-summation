! Copyright (C) 2020 Jonas A. Finkler
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module constants
    use precision
    implicit none

    real(dp), parameter :: ang2bohr = 1.8897259886_dp
    real(dp), parameter :: ev2ha    = 0.03674932539796232_dp
    real(dp), parameter :: pascal2atomic = 1._dp / 2.9421015697e13_dp ! hartree / bohr**3
    real(dp), parameter :: kBoltzmann = 3.166811965e-6_dp! atomic units (Ha/K)

    ! maximum element index
    integer, parameter :: maxElemNum = 118

    real(dp), parameter :: PI = 4 * atan (1.0_dp)

    character(len=2), parameter :: elemSymbols(118) = & ! from https://www.plaintextlist.com/science/list_of_chemical_elements_(symbols)
            [' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', ' S', &
             'Cl', 'Ar', ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', &
             'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
             'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', &
             'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
             'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', ' U', 'Np', 'Pu', 'Am', 'Cm', &
             'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', &
             'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

    ! these are the average weights for the isotope distribution as found on earth
    ! maybe the weight of the most common isotope would make more sense?
    real(dp), parameter :: atomicMasses(118) = & ! from https://www.qmul.ac.uk/sbcs/iupac/AtWt/index.html#02 (2019)
            [   1.008_dp,    4.002_dp,     6.94_dp,    9.012_dp,    10.81_dp,   12.011_dp,   14.007_dp,   15.999_dp, &
               18.998_dp,  20.1797_dp,   22.989_dp,   24.305_dp,   26.981_dp,   28.085_dp,   30.973_dp,    32.06_dp, &
                35.45_dp,   39.948_dp,  39.0983_dp,   40.078_dp,   44.955_dp,   47.867_dp,  50.9415_dp,  51.9961_dp, &
               54.938_dp,   55.845_dp,   58.933_dp,  58.6934_dp,   63.546_dp,    65.38_dp,   69.723_dp,   72.630_dp, &
               74.921_dp,   78.971_dp,   79.904_dp,   83.798_dp,  85.4678_dp,    87.62_dp,   88.905_dp,   91.224_dp, &
               92.906_dp,    95.95_dp,      97._dp,   101.07_dp,  102.905_dp,   106.42_dp, 107.8682_dp,  112.414_dp, &
              114.818_dp,  118.710_dp,  121.760_dp,   127.60_dp,  126.904_dp,  131.293_dp,  132.905_dp,  137.327_dp, &
              138.905_dp,  140.116_dp,  140.907_dp,  144.242_dp,     145._dp,   150.36_dp,  151.964_dp,   157.25_dp, &
              158.925_dp,  162.500_dp,  164.930_dp,  167.259_dp,  168.934_dp,  173.045_dp, 174.9668_dp,  178.486_dp, &
              180.947_dp,   183.84_dp,  186.207_dp,   190.23_dp,  192.217_dp,  195.084_dp,  196.966_dp,  200.592_dp, &
               204.38_dp,    207.2_dp,  208.980_dp,     209._dp,     210._dp,     222._dp,     223._dp,     226._dp, &
                 227._dp, 232.0377_dp,  231.035_dp,  238.028_dp,     237._dp,     244._dp,     243._dp,     247._dp, &
                 247._dp,     251._dp,     252._dp,     257._dp,     258._dp,     259._dp,     262._dp,     267._dp, &
                 270._dp,     269._dp,     270._dp,     270._dp,     278._dp,     281._dp,     281._dp,     285._dp, &
                 286._dp,     289._dp,     289._dp,     293._dp,     293._dp,     294._dp ]


contains

    elemental function getAtomicMass(n) result(R)
        integer, intent(in) :: n
        real(dp) :: R

        R = atomicMasses(n)

    end function

    function elemNumToSym(n) result(sym)
        integer, intent(in) :: n
        character(len=2) :: sym
        if (n>118) stop "no data for elements > 118 available"
        sym = elemSymbols(n)
    end function elemNumToSym

    function elemSymToNum(sym) result(n)
        character(len=*), intent(in) :: sym
        integer :: n
        integer :: i

        do i=1,118
            if (trim(adjustl(sym)) == trim(adjustl(elemSymbols(i)))) then
                n = i
                return
            end if
        end do
        print*, 'Error: Element ' // sym // ' not found.'
        stop
    end function elemSymToNum


end module constants


