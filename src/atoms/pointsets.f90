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


module pointsets
    use precision
    use linalg
    implicit none



contains

    !subroutine alignLatticeWithAxes(nat, lat, ats)
    !    integer, intent(in) :: nat
    !    real(dp), intent(inout) :: lat(3, 3), ats(3, nat)
    !    real(dp) :: newlat(3,3), R(3,3)
    !    integer :: i
    !    real(dp) :: rotdet

    !    ! calculate lat = R * newlat
    !    ! then newlat = R^t * lat
    !    call qrDecomposition(3, lat, R, newlat)
    !    rotdet = R(1,1) + R(2,2) + R(3,3)
    !    if (rotdet < 0._dp) then ! we have an improper rotation -> fix it
    !        R(:,:) = -1._dp * R(:,:)
    !        newlat(:,:) = -1._dp * newlat(:,:)
    !    end if

    !    R = transpose(R)

    !    do i=1, nat
    !        ats(:,i) = matMulVec3(R, ats(:,i))
    !    end do

    !    lat(:,:) = newlat(:,:)

    !end subroutine

    subroutine toRelativeCoordinates(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: ats(3, nat)
        real(dp), intent(in) :: lat(3,3)
        real(dp) :: invlat(3,3)
        integer :: i

        call inv3DM(lat, invlat)

        do i=1, nat
            ats(:,i) = mAtMulVec3(invlat, ats(:,i))
        end do

    end subroutine

    subroutine toAbsoluteCoordinates(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(inout) :: ats(3,nat)
        integer :: i

        do i=1, nat
            ats(:,i) = matMulVec3(lat, ats(:,i))
        end do

    end subroutine

    subroutine moveAtomsIntoCell(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(inout) :: ats(3, nat)
        integer :: i

        call toRelativeCoordinates(nat, lat, ats)

        do i=1, nat
            ats(1,i) = modulo(ats(1,i), 1._dp)
            ats(2,i) = modulo(ats(2,i), 1._dp)
            ats(3,i) = modulo(ats(3,i), 1._dp)
        end do

        call toAbsoluteCoordinates(nat, lat, ats)

    end subroutine

    ! returns the number of cells in the direction of each lattice vector that are needed to find all neighbors within cutoff
    subroutine getNcells(lattice, cutoff, n)
        implicit none
        real(dp), intent(in) :: lattice(3,3)
        real(dp), intent(in) :: cutoff
        integer, intent(out) :: n(3)
        real(dp) :: axb(3), axc(3), bxc(3)
        real(dp) :: proj(3)


        ! calculate the norm vectors of the 3 planes
        axb = cross(lattice(:,1), lattice(:,2))
        axb = axb / sqrt(sum(axb**2))
        axc = cross(lattice(:,1), lattice(:,3))
        axc = axc / sqrt(sum(axc**2))
        bxc = cross(lattice(:,2), lattice(:,3))
        bxc = bxc / sqrt(sum(bxc**2))

        proj(1) = dot(lattice(:,1), bxc(:))
        proj(2) = dot(lattice(:,2), axc(:))
        proj(3) = dot(lattice(:,3), axb(:))

        n(:) = ceiling(cutoff / abs(proj))

    end subroutine getNcells

    subroutine translate(nat, atoms,r)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: atoms(3, nat)
        real(dp), intent(in) :: r(3)
        integer :: i

        do i=1,nat
            atoms(:,i) = atoms(:,i) + r(:)
        end do
    end subroutine

    !points are gaussian distributed. (this is quite expensive (lots of math))
    subroutine randomPointSet(nat, ats)
        use random
        integer, intent(in) :: nat
        real(dp), intent(out) :: ats(3, nat)
        integer :: i, j

        do j=1,nat
            do i=1,3
                ats(i,j) = random_normal()
            end do
        end do

    end subroutine randomPointSet


end module

