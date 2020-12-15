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


module linalg
    use precision
    use constants

    implicit none

contains

    ! only for symmetric matrices so far
    !subroutine qrDecomposition(n, A, Q, R)
    !    integer, intent(in) :: n
    !    real(dp), intent(in) :: A(n,n)
    !    real(dp), intent(out) :: Q(n,n), R(n,n)

    !    real(dp) :: tau(n), workq(1)
    !    real(dp), allocatable :: work(:)
    !    integer :: lwork, info, i

    !    R = A
    !    call dgeqrf(n, n, R, n, tau, workq, -1, info)
    !    lwork = int(workq(1))
    !    allocate(work(int(lwork)))
    !    call dgeqrf(n, n, R, n, tau, work, lwork, info)
    !    if (info /= 0) stop 'QR decomposition failed!'
    !    deallocate(work)

    !    Q = 0._dp
    !    do i=1,n
    !        Q(i,i) = 1._dp
    !    end do

    !    call dormqr('L', 'N', n, n, n, R, n, tau, Q, n, workq, -1, info)
    !    lwork = int(workq(1))
    !    allocate(work(int(lwork)))
    !    call dormqr('L', 'N', n, n, n, R, n, tau, Q, n, work, lwork, info)
    !    if (info /= 0) stop 'QR decomposition failed! (2)'
    !    deallocate(work)

    !    ! lapack stores information about Q in the lower half of R. Removing it.
    !    do i=1,n-1
    !        R(i+1:n,i) = 0._dp
    !    end do

    !end subroutine QRdecomposition

    function vecMulVecT(n, a, b) result(c)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
        real(dp) :: c(n,n)
        integer :: i, j

        do i=1,n
            do j=1,n
                c(i,j) = a(i) * b(j)
            end do
        end do


    end function

    function matMulVec3(mat, vec) result(res)
        implicit none
        real(dp), dimension(3, 3), intent(in) :: mat
        real(dp), dimension(3), intent(in) :: vec
        real(dp), dimension(3) :: res

        res(1) = sum(mat(1, :) * vec(:))
        res(2) = sum(mat(2, :) * vec(:))
        res(3) = sum(mat(3, :) * vec(:))
    end function

    function matMulVecNM(n, m, mat, vec) result(res)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: mat(n,m), vec(m)
        real(dp) :: res(n)
        integer :: i

        do i=1,n
            res(i) = sum(mat(i,:) * vec(:))
        end do

    end function

    function cross(a, b) result(c)
        implicit none
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3) :: c

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function

    function det3D(M) result(det)
        real(dp), intent(in) :: M(3,3)
        real(dp) :: det

        det = M(1,1) * M(2,2) * M(3,3) &
            + M(1,2) * M(2,3) * M(3,1) &
            + M(1,3) * M(2,1) * M(3,2) &
            - M(1,1) * M(2,3) * M(3,2) &
            - M(1,2) * M(2,1) * M(3,3) &
            - M(1,3) * M(2,2) * M(3,1)

    end function det3D

    subroutine inv3D(a1, a2, a3, b1, b2, b3)
        implicit none

        real(dp), intent(in), dimension(3) :: a1, a2, a3
        real(dp), intent(out), dimension(3) :: b1, b2, b3
        real(dp) :: det

        det = a1(1) * a2(2) * a3(3) &
            + a1(2) * a2(3) * a3(1) &
            + a1(3) * a2(1) * a3(2) &
            - a1(1) * a2(3) * a3(2) &
            - a1(2) * a2(1) * a3(3) &
            - a1(3) * a2(2) * a3(1)

        b1 = [a2(2) * a3(3) - a2(3) * a3(2), a1(3) * a3(2) - a1(2) * a3(3), a1(2) * a2(3) - a1(3) * a2(2)]
        b2 = [a2(3) * a3(1) - a2(1) * a3(3), a1(1) * a3(3) - a1(3) * a3(1), a1(3) * a2(1) - a1(1) * a2(3)]
        b3 = [a2(1) * a3(2) - a2(2) * a3(1), a1(2) * a3(1) - a1(1) * a3(2), a1(1) * a2(2) - a1(2) * a2(1)]

        b1 = b1 / det
        b2 = b2 / det
        b3 = b3 / det
    end subroutine

    subroutine inv3DM(A, B)
        implicit none
        real(dp), intent(in), dimension(3, 3) :: A
        real(dp), intent(out), dimension(3, 3) :: B
        real(dp), dimension(3) :: a1, a2, a3
        real(dp), dimension(3) :: b1, b2, b3
        real(dp) :: det

        a1 = A(1, :)
        a2 = A(2, :)
        a3 = A(3, :)

        det = det3D(A)

        b1 = [a2(2) * a3(3) - a2(3) * a3(2), a1(3) * a3(2) - a1(2) * a3(3), a1(2) * a2(3) - a1(3) * a2(2)]
        b2 = [a2(3) * a3(1) - a2(1) * a3(3), a1(1) * a3(3) - a1(3) * a3(1), a1(3) * a2(1) - a1(1) * a2(3)]
        b3 = [a2(1) * a3(2) - a2(2) * a3(1), a1(2) * a3(1) - a1(1) * a3(2), a1(1) * a2(2) - a1(2) * a2(1)]

        b1 = b1 / det
        b2 = b2 / det
        b3 = b3 / det

        B(1, :) = b1
        B(2, :) = b2
        B(3, :) = b3
    end subroutine

    function dot(a,b)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: dot
        if (size(a) /= size(b)) then
            print*, 'dot product uf array with not same length'
            stop
        end if
        dot = sum(a(:) *  b(:))
    end function

    subroutine recLattice(lat, reclat)
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(out) :: reclat(3,3)

        call inv3DM(transpose(lat), reclat)
        reclat(:,:) = reclat(:,:) * (2._dp * PI )

    end subroutine recLattice

    function vecAngle(a, b) result(t)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: t

        t = acos(sum(a*b) / sqrt(sum(a**2) * sum(b**2)))

    end function vecAngle

end module linalg