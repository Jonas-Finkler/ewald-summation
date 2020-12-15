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


! this is the place to put small utility functions
module util
    implicit none

contains

    subroutine getDistinctElements(n, A, nel, els, count)
        integer, intent(in) :: n
        integer, intent(in) :: A(n)
        integer, intent(out) :: nel
        integer, intent(out), allocatable :: els(:)
        integer, intent(out), optional, allocatable :: count(:)

        integer :: i, j
        integer :: tmp(n), tmpCount(n)
        logical :: found

        tmp(:) = -1
        tmpCount(:) = 0
        nel = 0
        do i=1,n
            found = .false.
            do j=1,nel
                if (A(i) == tmp(j)) then
                    found = .true.
                    tmpCount(j) = tmpCount(j) + 1
                    exit
                end if
            end do
            if (.not. found) then
                nel = nel + 1
                tmp(nel) = A(i)
                tmpCount(nel) = 1
            end if
        end do

        if (allocated(els)) deallocate(els)
        allocate(els(nel))
        els(:) = tmp(:nel)

        if (present(count)) then
            if(allocated(count)) deallocate(count)
            allocate(count(nel))
            count(:) = tmpCount(:nel)
        end if

    end subroutine getDistinctElements

    function arrayCompare(a, b)
        logical :: arrayCompare
        integer, intent(in) :: a(:), b(:)
        integer :: i
        arrayCompare = .false.
        if(size(a) /= size(b)) return
        do i = 1, size(a)
            if(a(i)/=b(i)) return
        end do
        arrayCompare = .true.
    end function arrayCompare

    subroutine assert(condition, errorMessage)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: errorMessage

        if(.not. condition) then
            write(*,*) 'assertion failed: ', errorMessage
            stop 
        end if

    end subroutine assert

end module util

