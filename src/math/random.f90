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



module random
    use precision
    implicit none

contains

    ! knuth shuffle
    subroutine shuffle(n, a)
        integer, intent(in) :: n
        integer, intent(inout) :: a(n)
        integer :: i, randpos, temp
        real(dp) :: r

        do i = n, 2, -1
            call random_number(r)
            randpos = int(r * i) + 1
            temp = a(randpos)
            a(randpos) = a(i)
            a(i) = temp
        end do

    end subroutine shuffle

    subroutine random_normal_1D(A)
        real(dp), intent(out) :: A(:)
        integer :: i
        integer :: n

        n = size(A)
        do i=1,n
            A(i) = random_normal()
        end do
    end subroutine

    subroutine random_normal_2D(A)
        real(dp), intent(out) :: A(:,:)
        integer :: i, j
        integer :: n(2)

        n = shape(A)
        do j=1,n(2)
            do i=1,n(1)
                A(i,j) = random_normal()
            end do
        end do
    end subroutine

    subroutine random_normal_3D(A)
        real(dp), intent(out) :: A(:,:,:)
        integer :: i, j, k
        integer :: n(3)

        n = shape(A)
        do k=1,n(3)
            do j=1,n(2)
                do i=1,n(1)
                    A(i,j,k) = random_normal()
                end do
            end do
        end do
    end subroutine

    ! uses the Box-Müller algorithm to generate normal distributed random numbers
    ! threadsafe -> can be used in omp parallelized sections
    function random_normal() result (rnor)
        real(dp), parameter :: PI = 4._dp * datan(1._dp)
        real(dp) :: rnor
        real(dp) :: u(2)
        real(dp), save :: cache
        logical, save :: has_cache = .false.
        ! each thread has its own cache
        !$omp threadprivate(cache, has_cache)

        if (has_cache) then
            has_cache = .false.
            rnor = cache
        else
            call random_number(u)
            rnor  = sqrt(-2._dp * log(u(1))) * sin(2._dp * PI * u(2))
            cache = sqrt(-2._dp * log(u(1))) * cos(2._dp * PI * u(2))
            has_cache = .true.
        end if
    end function random_normal

        ! xorshift64*
    ! took implementation of the algorithm in the Spire package as reference (MIT license)
    ! https://github.com/typelevel/spire/blob/master/extras/src/main/scala/spire/random/rng/XorShift64Star.scala
    ! has a period of 2**64-1
    ! translated to fortran by JAF
    function xorshift64star(state) result(rnd)
        implicit none
        real*8 :: rnd
        integer*8, intent(inout) :: state
        ! integer(16) :: x
        ! integer*8 :: xx
        integer*8 :: b(4)
        integer*8, parameter :: f(4) = [56605, 20332, 62609, 9541] ! 2685821657736338717 in base 2**16
        integer*8 :: m(4)

        state = ieor(state, ishft(state, -12))
        state = ieor(state, ishft(state,  25))
        state = ieor(state, ishft(state, -27))

        ! fortran has no unsigned ints.
        ! below would be how we convert signed to unsigned (twos complement)
        ! if (x < 0) then
        !   x = 2_16**64 + state
        ! else
        !   x = state
        ! end if
        ! doing the modulo is equivalent however
        ! BUT: no 128 bit int support from intel compiler
        ! we therefore use some trickery by transforming the number into base 2**16
        ! the version below would work with gnu compiler
        ! x = modulo(state * 2685821657736338717_16, 2_16**64)
        ! rnd = real(x, 8) / 2_16**64
        ! write(*,*) rnd
        ! this is the version for intel
        call toBase16(state, b)
        call multBase16(b, f, m)
        ! call fromBase16(m, xx)
        ! to return a double we divide by 2^64
        rnd = m(4) / 2.d0**16 + m(3) / 2.d0**32 + m(2) / 2.d0**48 + m(1) / 2.d0**64
        ! write(*,*) rnd

    end function

    ! becaus we cannot use 128 bit ints in ifort we represent our number using 4 ints in base 2**16
    ! x = sum_i b(i) * 2**(16*(i-1))
    ! treats x as unsigned int
    subroutine toBase16(x, b)
        implicit none
        integer*8, intent(in) :: x
        integer*8, intent(out) :: b(4)
        integer :: i
        integer*8 :: t

        b(:) = 0
        if (x<0) then
            t = (9223372036854775807_8 + x) + 1
            b(4) = 2_8**15
        else
          t = x
        end if

        do i=1,4
            b(i) = b(i) + modulo(t, 2_8**16)
            t = t / 2_8**16
        end do

    end subroutine

    ! also treats x as unsigned int
    subroutine fromBase16(b, x)
        implicit none
        integer*8, intent(in) :: b(4)
        integer*8, intent(out) :: x
        integer :: i

        x = 0
        do i=1,3
            x = x + b(i) * 2_8**(16*(i-1))
        end do
        x = x + modulo(b(4), 2_8**16) * 2_8**(16*3)
    end subroutine

    ! does multiplication mod 2**64 of two numbers in the base of 2**16
    subroutine multBase16(a, b, c)
        implicit none
        integer*8, intent(in) :: a(4), b(4)
        integer*8, intent(out) :: c(4)

        integer :: i,j

        c(:) = 0
        do i=1,4
            do j=1,5-i
                c(i+j-1) = c(i+j-1) + a(i) * b(j)
            end do
        end do
        do i=1,3
            c(i+1) = c(i+1) + c(i) / 2_8**16
            c(i) = modulo(c(i), 2_8**16)
        end do
        c(4) = modulo(c(4), 2_8**16)

    end subroutine


end module random
