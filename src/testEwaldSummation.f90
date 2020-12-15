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


program testEwaldSummation
    use precision
    use ewaldSummation
    use linalg
    use pointsets
    use parameters
    implicit none

    integer, parameter :: nat = 10
    real(dp), parameter :: dx = 1.e-5_dp
    real(dp), parameter :: refMadelung = -1.7475645946331821906362120355443_dp
    logical, parameter :: pointCharges = .false.

    real(dp) :: sigma(nat)
    real(dp) :: ats(3,nat)
    real(dp) :: q(nat)
    real(dp) :: lat(3,3), dlat(3,3)
    real(dp) :: e1, f1(3,nat)
    real(dp) :: e2, f2(3,nat)
    real(dp) :: eps(3,3), stressfd(3,3), stressan(3,3)
    real(dp) :: fracc(3,nat), tm(3,3), unit(3,3), tmplat(3,3)
    real(dp) :: disp(3,nat)
    integer :: i, j
    real(dp) :: A(nat, nat), dAdxyzq(nat, 3, nat)
    real(dp) :: zeros(nat)

    real(dp) :: atsMadelung(3,8), eMadelung, fMadelung(3,8), qMadelung(8), sigmaMadelung(8)


    print*, 'Welcome to the Ewald summation test routine.'
    print*, 'This code was written by Jonas A. Finkler and was published under the GNU General Public License v3.'
    print*, 'The source code can be found here: https://github.com/Jonas-Finkler/ewald-summation'
    print*, ''

    print*, 'The accuracy was set to: ', ewaldSummationPrecision
    print*, ''


    ! unit matrix
    unit = 0._dp; unit(1,1) = 1._dp; unit(2,2) = 1._dp; unit(3,3) = 1._dp
    zeros(:) = 0._dp

    ! Lets test if we get the correct result for the Madelung Constant
    lat(:,:) = 0._dp
    lat(1,1) = 2._dp
    lat(2,2) = 2._dp
    lat(3,3) = 2._dp
    ! Na
    atsMadelung(:,1) = [0._dp, 0._dp, 0._dp]
    atsMadelung(:,2) = [1._dp, 1._dp, 0._dp]
    atsMadelung(:,3) = [1._dp, 0._dp, 1._dp]
    atsMadelung(:,4) = [0._dp, 1._dp, 1._dp]
    ! Cl
    atsMadelung(:,5) = [1._dp, 1._dp, 1._dp]
    atsMadelung(:,6) = [1._dp, 0._dp, 0._dp]
    atsMadelung(:,7) = [0._dp, 1._dp, 0._dp]
    atsMadelung(:,8) = [0._dp, 0._dp, 1._dp]

    ! Na
    qMadelung(:4) = 1._dp
    ! Cl
    qMadelung(5:) = -1._dp

    ! Use point charges
    sigmaMadelung = -1._dp

    call ewaldEnergy(8, atsMadelung, lat, qMadelung, sigmaMadelung, eMadelung, fMadelung, stressan)
    eMadelung = eMadelung / 4._dp ! 4 atoms of each type

    print*, 'Testing energy using Madelung constant'
    print*, 'Reference value of the NaCl Madelung Constant:'
    print '(F25.20)', refMadelung
    print*, 'Result obtained using Ewald summation:'
    print '(F25.20)', eMadelung
    print*, 'Difference (ref-ewald):'
    print '(E25.18)', refMadelung - eMadelung
    print*, 'Total force (should be zero):'
    print*, sqrt(sum(fMadelung**2))
    print*, ''



    print*, 'Now the force and stress tensor will be tested using finite differences.'
    print*, 'Number of atoms: ', nat
    print*, 'Using point charges: ', pointCharges

    ! generate some random lattice
    lat(:,1) = [1.0_dp, 0.0_dp, 0.0_dp]
    lat(:,2) = [0.0_dp, 1.0_dp, 0.0_dp]
    lat(:,3) = [0.0_dp, 0.0_dp, 1.0_dp]
    call random_number(dlat)
    lat = lat + dlat

    ! random atomic positions
    call random_number(ats)
    do i=1,nat
        ats(:,i) = matMulVec3(lat, ats(:,i))
    end do

    ! random charges that add up to 0
    call random_number(q)
    q = q - sum(q) / real(nat, dp)

    if (pointCharges) then
    ! random Gaussian charges of std. dev. sigma
        call random_number(sigma)
        sigma(:) = sigma(:) * 0.03_dp + 0.07_dp
    ! point charges
    else
        sigma(:) = -1._dp
    end if

    ! forces
    ! take some random displacement
    call random_number(disp)
    disp = disp - 0.5_dp
    disp = disp * dx
    ! calculate change of energy along displacememnt
    call ewaldEnergy(nat, ats, lat, q, sigma, e1, f1)
    fracc = ats + disp
    call ewaldEnergy(nat, fracc, lat, q, sigma, e2, f2)

    ! compare to result that we get using the force from the Ewald sum
    print*, 'Energy change after a small displacement (A):'
    print*, e2 - e1
    print*, 'Result obtained from dE=-sum(F*dx) (B):'
    print*, -1._dp * sum(disp * (f1+f2) / 2._dp)
    print*, 'Difference (A-B)'
    print*, (e2 - e1) - (-1._dp * sum(disp * (f1+f2) / 2._dp))
    print*, 'Ratio (A/B)'
    print*, (e2 - e1) / (-1._dp * sum(disp * (f1+f2) / 2._dp))
    print*, ''



    ! Calculate the anayltic stress from the ewald sum
    call ewaldEnergy(nat, ats, lat, q, sigma, e1, f1, stressan)

    ! Calculate finite difference stress
    do i=1,3
        do j=1,3
            fracc = ats
            call toRelativeCoordinates(nat, lat, fracc)
            eps = 0._dp
            eps(i,j) = dx
            tm = unit + eps
            tmplat = matmul(tm, lat)
            call toAbsoluteCoordinates(nat, tmplat, fracc)
            call ewaldEnergy(nat, fracc, tmplat, q, sigma, e2, f2)
            stressfd(i,j) = (e2 - e1) / dx
        end do
    end do
    stressfd = -1._dp / det3D(lat) * stressfd

    print*, 'Analytic stress tensor obtained from Ewald summation (A):'
    do i=1,3
        print*, stressan(i,:)
    end do
    print*, 'Result obtained using finite differences (B):'
    do i=1,3
        print*, stressfd(i,:)
    end do

    print*, 'Difference (A-B)'
    do i=1,3
        print*, stressan(i,:) - stressfd(i,:)
    end do

    print*, 'Ration (A/B)'
    do i=1,3
        print*, stressan(i,:) / stressfd(i,:)
    end do

    print*, ''

    ! We can also calculate the energy using the A Matrix
    ! But be aware that this is less efficient and not as optimized
    call ewaldEnergy(nat, ats, lat, q, sigma, e1, f1)
    call eemMatrixEwald(nat, ats, lat, sigma, zeros, A)
    e2 = 0.5_dp * sum(q * matmul(A, q))

    print*, 'The energy can also be calculate using the Matrix '
    print*, 'Energy from A Matrix:'
    print*, e2
    print*, 'Difference to direct result:'
    print*, e2-e1
    print*, 'Ratio to direct result:'
    print*, e2 / e1
    print*, ''

    ! We can also get the forces
    call eemdAdxyzTimesQEwald(nat, ats, lat, sigma, q, dAdxyzQ)
    f2(:,:) = 0._dp
    do i=1,nat
        f2(:,:) = f2(:,:) + (-0.5_dp) * q(i) * dAdxyzq(i,:,:)
    end do
    print*, 'The forces can also be calculated using the matrix dAdxyzq'
    print*, 'Absolute force error to direct result:'
    print*, sum((f2-f1)**2)
    print*, ''

end program testEwaldSummation