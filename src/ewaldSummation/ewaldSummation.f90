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



! stress: https://wanglab.hosted.uark.edu/DLPOLY2/node65.html

! IMPORTANT: if the first element of sigma is set to something <= 0, point charges are used
! in this case self interaction of the charges is ommitted

! lattice convention: lat(:,i) is the i-th lattice vector

module ewaldSummation
    use precision
    use constants
    use pointsets
    use linalg
    use parameters

    implicit none
    ! In RuNNer: ewaldalpha = 1 / (2 * eta)
    private
    public :: ewaldEnergy,  eemMatrixEwald, eemdAdxyzTimesQEwald, dEdLatFromStress


contains


    ! https://arxiv.org/pdf/1904.08875.pdf
    ! https://aip.scitation.org/doi/pdf/10.1063/1.5053968
    ! Be careful: in the first paper alpha != alpha_{RuNNer}
    ! alpha = 1 / (sqrt(2) eta)
    function getOptimalEta(nat, lat) result(eta)
        integer, intent(in) :: nat
        real(dp), intent(in) :: lat(3,3)
        real(dp) :: eta

        eta = 1._dp / sqrt(2._dp * PI) * (det3D(lat)**2 / nat)**(1._dp / 6._dp)

    end function getOptimalEta

    ! When the matrix is built, each reciprocal term takes O(nat**2) compuations
    ! Therefore reciprocal is not cheaper than real anymore and the optimal choice of eta is different (it is independant of nat)
    function getOptimalEtaMatrix(lat) result(eta)
        real(dp), intent(in) :: lat(3,3)
        real(dp) :: eta

        eta = 1._dp / sqrt(2._dp * PI) * det3D(lat)**(1._dp / 3._dp)

    end function getOptimalEtaMatrix

    function getOptimalCutoffReal(eta, A) result(r)
        real(dp), intent(in) :: eta
        real(dp), intent(in) :: A
        real(dp) :: r

        r = sqrt(2._dp) * eta * sqrt(-log(A))

    end function getOptimalCutoffReal

    function getOptimalCutoffRecip(eta, A) result(r)
        real(dp), intent(in) :: eta
        real(dp), intent(in) :: A
        real(dp) :: r

        r = sqrt(2._dp) / eta * sqrt(-log(A))

    end function getOptimalCutoffRecip

    subroutine dEdlatFromStress(lat, stress, dEdlat)
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(in) :: stress(3,3)
        real(dp), intent(out) :: dEdlat(3,3)
        real(dp) :: invlat(3,3)

        call inv3DM(lat, invlat)
        dEdlat = -1._dp * det3D(lat) * matmul(stress, transpose(invlat))

    end subroutine

    ! gaussian charges (and also ewald gaussians) are of the form 1/(2*pi*sigma**2)**(3/2) * exp(-r**2/(2*sigma**2))
    ! -> sigma = sqrt(2) * alpha_{cent}
    subroutine ewaldEnergy(nat, ats, lat, q, sigma, e, f, stress)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), q(nat), sigma(nat)
        real(dp), intent(out) :: e, f(3,nat)
        real(dp), intent(out), optional :: stress(3,3)
        real(dp) :: ereal, erecip, eself
        real(dp ) :: stressreal(3,3), stressrecip(3,3)
        real(dp) :: freal(3,nat), frecip(3,nat)
        real(dp) :: eta

       ! real(dp) :: et,st

        eta = getOptimalEta(nat, lat)
        eta = max(eta, maxval(sigma))

       ! call cpu_time(st)
        call ewaldSumReal(nat, ats, lat, q, sigma, eta, ereal, freal, stressreal)
       ! call cpu_time(et)
       ! print*, 'time real', et-st
       ! call cpu_time(st)
        call ewaldSumRecip(nat, ats, lat, q, eta, erecip, frecip, stressrecip)
       ! call cpu_time(et)
       ! print*, 'time recip', et-st
        call ewaldSelf(nat, q, eta, eself)

        e = ereal + erecip + eself
        f = freal + frecip

        if (present(stress)) then
            stress = (stressreal + stressrecip) / det3D(lat)
        end if

    end subroutine ewaldEnergy

    subroutine ewaldSumReal(nat, ats, lat, q, sigma, eta, e, f, stress)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), q(nat), sigma(nat)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: e, f(3,nat), stress(3,3)

        real(dp) :: d(3), r, dlat(3), inv2eta, invsqrt2eta, interf, gamma, d2, cutoff, ftmp(3)
        integer :: i, j, k, iat, jat, n(3), sj, sk, jats

        invsqrt2eta = 1._dp / (sqrt(2._dp) * eta)
        inv2eta = 1._dp / (2._dp * eta)
        e = 0._dp
        f(:,:) = 0._dp
        stress(:,:) = 0._dp


        cutoff = getOptimalCutoffReal(eta, ewaldSummationPrecision)
        call getNcells(lat, cutoff, n)

        !$omp parallel do private(i, sj, j, sk, k, dlat, jats, d, d2, r, interf, gamma, ftmp)&
        !$omp reduction(+:e,f,stress)
        do iat=1,nat
            do i=0,n(1)
                sj = -n(2)
                if (i==0) sj = 0
                do j=sj,n(2)
                    sk = -n(3)
                    if (i==0.and.j==0) sk = 0
                    do k=sk,n(3)
                        dlat(:) = i * lat(:,1) + j * lat(:,2) + k * lat(:,3)
                        jats = 1
                        if (i==0 .and. j==0 .and. k==0) jats=iat+1
                        do jat=jats,nat
                            if (i/=0 .or. j/=0 .or. k/=0 .or. iat/=jat) then
                                d(:) = ats(:,iat) - ats(:,jat) + dlat(:)
                                d2 = sum(d**2)
                                if (d2 > cutoff**2) cycle
                                r = sqrt(d2)
                                interf = erfc(r * invsqrt2eta)
                                if (sigma(1) > 0._dp) then
                                    gamma = sqrt(sigma(iat)**2 + sigma(jat)**2)
                                    !interf = erf(r / (sqrt(2._dp) * gamma)) - erf(r * inv2eta)
                                    interf = interf - erfc(r / (sqrt(2._dp) * gamma))
                                end if

                                e = e + 2._dp * q(iat) * q(jat) * interf / r

                                if (sigma(1) <= 0._dp) then
                                    ftmp(:) = &
                                      + q(iat) * q(jat) * d(:) / r**3 &
                                      * (2._dp * r * invsqrt2eta * exp(-invsqrt2eta**2 * r**2) / sqrt(PI) &
                                      + interf)
                                else
                                    ftmp(:) = &
                                      + q(iat) * q(jat) * d(:) / r**3 &
                                      * (2._dp * r * invsqrt2eta * exp(-invsqrt2eta**2 * r**2) / sqrt(PI) &
                                      - 2._dp * r / (sqrt(2._dp) * gamma) * exp(-1._dp / (2._dp * gamma**2) * r**2) / sqrt(PI) &
                                      + interf)
                                end if
                                f(:,iat) = f(:,iat) + ftmp(:)
                                f(:,jat) = f(:,jat) - ftmp(:)
                                stress(:,:) = stress(:,:) + vecMulVecT(3, d, ftmp)
                            end if
                        end do
                    end do
                end do
            end do
            if (sigma(1) > 0._dp) then ! self interaction
                e = e + q(iat)**2 / (sigma(iat) * 2._dp * sqrt(PI)) * 2._dp
            end if
        end do
        !$omp end parallel do
        e = e / 2._dp

    end subroutine ewaldSumReal

    subroutine ewaldSumRecip(nat, ats, lat, q, eta, e, f, stress)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), q(nat)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: e, f(3,nat), stress(3,3)

        real(dp) :: reclat(3,3), dlat(3), V, r, Sreal, Simag, kr, dSreal, dSimag, factor, cutoff
        integer :: i, j, k, iat, n(3), sj, sk, nklats, ik
        integer, allocatable :: klats(:,:)
        real(dp), parameter :: unit(3,3) = reshape([1._dp, 0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 1._dp], [3,3])

        call recLattice(lat, reclat)
        V = det3D(lat)
        e = 0._dp
        f = 0._dp
        stress = 0._dp

        cutoff = getOptimalCutoffRecip(eta, ewaldSummationPrecision)
        call getNcells(reclat, cutoff, n)


        allocate(klats(3,(n(1)+1)*(n(2)*2+1)*(n(3)*2+1)))
        nklats = 0


        ! build a list of all kpoints so that we can parallelize the loop over them
        do i=0,n(1)
            sj = -n(2)
            if (i==0) sj = 0
            do j=sj,n(2)
                sk = -n(3)
                if (i==0.and.j==0) sk = 0
                do k=sk,n(3)
                    nklats = nklats + 1
                    if (nklats>(n(1)+1)*(n(2)*2+1)*(n(3)*2+1)) stop 'too many k points'
                    klats(:,nklats) = [i,j,k]
                end do
            end do
        end do


        !$omp parallel do default(shared) private(i,j,k,dlat,r,Sreal,Simag,kr,factor,iat,dSreal,dSimag) &
        !$omp reduction(+:e,f,stress)
        do ik=1,nklats
            i = klats(1,ik)
            j = klats(2,ik)
            k = klats(3,ik)
            if (i/=0 .or. j/=0 .or. k/=0) then
                dlat(:) = i * reclat(:,1) + j * reclat(:,2) + k * reclat(:,3)
                r = sum(dlat**2) ! r = k**2
                if (r > cutoff**2) cycle
                factor = exp(-eta**2 * r / 2._dp) / r
                Sreal = 0._dp
                Simag = 0._dp
                do iat=1,nat
                    kr = sum(dlat(:) * ats(:,iat))
                    Sreal = Sreal + q(iat) * cos(kr)
                    Simag = Simag + q(iat) * sin(kr)
                end do
                e = e + 2._dp * factor * (Sreal**2 + Simag**2)
                stress(:,:) = stress(:,:) + (unit - vecMulVecT(3, dlat, dlat) *(eta**2 + 2._dp / r)) &
                        * 2._dp * factor * (Sreal**2 + Simag**2)
                do iat=1,nat
                    kr = sum(dlat(:) * ats(:,iat))
                    dSreal = q(iat) * (-sin(kr))
                    dSimag = q(iat) * cos(kr)

                    f(:,iat) = f(:,iat) - 2._dp * dlat(:) * factor &
                            * 2._dp * (Sreal * dSreal + Simag * dSimag)

                end do
            end if
        end do
        !$omp end parallel do

        e =          e  / V * 4._dp * PI  / 2._dp
        f(:,:) = f(:,:) / V * 4._dp * PI  / 2._dp
        stress(:,:) = stress(:,:) / V * 4._dp * PI  / 2._dp

    end subroutine ewaldSumRecip

    subroutine ewaldSelf(nat, q, eta, e)
        integer, intent(in) :: nat
        real(dp), intent(in) :: q(nat)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: e
        e = -sum(q**2) / (sqrt(2._dp * PI) * eta)
    end subroutine ewaldSelf

    ! #######################################################################################################

    ! Calculation of Matrix A, sucht that E = 1/2 * q^t * A * q

    ! the periodic A matrix is not unique (depends on eta)
    ! however, because sum(q)==0, A*q is independent of eta
    ! todo: These subroutines are not optimized and parallelized yet
    subroutine eemMatrixEwald(nat, ats, lat, sigma, hardness, A)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), sigma(nat), hardness(nat)
        real(dp), intent(out) :: A(nat, nat)
        real(dp) :: Atmp(nat, nat), eta
        !        real(dp) :: et, st

        eta = getOptimalEtaMatrix(lat)
        eta = max(eta, maxval(sigma))

        A(:,:) = 0._dp
        !        call cpu_time(st)
        call AmatrixReal(nat, ats, lat, sigma, hardness, eta, Atmp)
        !        call cpu_time(et)
        !        write(*,*) 'real', et-st
        A(:,:) = A(:,:) + Atmp(:,:)
        !        call cpu_time(st)
        call AmatrixRecip(nat, ats, lat, eta, Atmp)
        !        call cpu_time(et)
        !        write(*,*) 'recip', et-st
        A(:,:) = A(:,:) + Atmp(:,:)
        !        call cpu_time(st)
        call AmatrixSelf(nat, eta, Atmp)
        !        call cpu_time(et)
        !        write(*,*) 'self', et-st
        A(:,:) = A(:,:) + Atmp(:,:)

    end subroutine eemMatrixEwald

    subroutine AmatrixReal(nat, ats, lat, sigma, hardness, eta, A)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), sigma(nat), hardness(nat)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: A(nat, nat)

        real(dp) ::  r, d(3), dlat(3), inv2eta, invsqrt2eta, gamma, interf, d2, cutoff
        integer :: i, j, k, iat, jat,n(3)

        invsqrt2eta = 1._dp / (sqrt(2._dp) * eta)
        inv2eta = 1._dp / (2._dp * eta)

        A = 0._dp

        cutoff = getOptimalCutoffReal(eta, ewaldSummationPrecision)
        call getNcells(lat, cutoff, n)

        do iat=1,nat
            do i=-n(1),n(1)
                do j=-n(2),n(2)
                    do k=-n(3),n(3)
                        dlat(:) = i * lat(:,1) + j * lat(:,2) + k * lat(:,3)
                        do jat=iat,nat
                            if (i/=0 .or. j/=0 .or. k/=0 .or. iat/=jat) then
                                d(:) = ats(:,iat) - ats(:,jat) + dlat(:)
                                d2 = sum(d**2)
                                if (d2 > cutoff**2) cycle
                                r = sqrt(d2)
                                interf = erfc(r * invsqrt2eta)
                                if (sigma(1) > 0._dp) then
                                    gamma = sqrt(sigma(iat)**2 + sigma(jat)**2)
                                    interf = interf - erfc(r / (sqrt(2._dp) * gamma))
                                end if
                                A(iat, jat) = A(iat, jat) + interf / r
                            end if
                        end do
                    end do
                end do
            end do
            if (sigma(1) > 0._dp) then ! self interaction
                A(iat, iat) = A(iat, iat) + 1._dp / (sigma(iat) * sqrt(PI)) + hardness(iat)
            end if
        end do

        do iat=1,nat
            do jat=iat+1,nat
                A(jat, iat) = A(iat,jat)
            end do
        end do

    end subroutine AmatrixReal

    ! please not that this routine has a different scaling that than the one ewaldEnergyRecip,
    ! as each reciprocal lattice point requires O(nat**2) instead of O(nat) computations
    ! hence the different optimal ewald eta for the matrix
    subroutine AmatrixRecip(nat, ats, lat, eta, A)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: A(nat, nat)
        real(dp) :: reclat(3,3), dlat(3), V, r, kri, krj, cutoff
        real(dp) :: factor
        integer :: i, j, k, iat, jat, n(3)

        call recLattice(lat, reclat)
        V = det3D(lat)

        A = 0._dp

        cutoff = getOptimalCutoffRecip(eta, ewaldSummationPrecision)
        call getNcells(reclat, cutoff, n)

        do i=-n(1),n(1)
            do j=-n(2),n(2)
                do k=-n(3),n(3)
                    if (i/=0 .or. j/=0 .or. k/=0) then
                        dlat(:) = i * reclat(:,1) + j * reclat(:,2) + k * reclat(:,3)
                        r = sum(dlat**2) ! r = k**2
                        if (r > cutoff**2) cycle
                        factor = exp(-eta**2 * r / 2._dp) / r
                        do iat=1,nat
                            kri = sum(dlat(:) * ats(:,iat))
                            do jat=iat,nat
                                krj = sum(dlat(:) * ats(:,jat))
                                A(iat,jat) = A(iat,jat) + factor * (cos(kri) * cos(krj) + sin(kri) * sin(krj))
                            end do
                        end do
                    end if
                end do
            end do
        end do

        do iat=1,nat
            do jat=iat+1,nat
                A(jat, iat) = A(iat,jat)
            end do
        end do

        A = A / V * 4._dp * PI

    end subroutine AmatrixRecip

    subroutine AmatrixSelf(nat, eta, A)
        integer, intent(in) :: nat
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: A(nat, nat)
        integer :: i

        A(:,:) = 0._dp
        do i=1,nat
            A(i,i) = -2._dp / (sqrt(2._dp * PI) * eta)
        end do

    end subroutine AmatrixSelf

    subroutine eemdAdxyzTimesQEwald(nat, ats, lat, sigma, q, dAdxyzQ)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), sigma(nat), q(nat)
        real(dp), intent(out) :: dAdxyzQ(nat, 3, nat)
        real(dp) :: tmp(nat, 3, nat), eta

        eta = getOptimalEtaMatrix(lat)
        eta = max(eta, maxval(sigma))

        call dAdxyzQReal(nat, ats, lat, sigma, q, eta, dAdxyzQ)
        call dAdxyzQRecip(nat, ats, lat, q, eta, tmp)
        dAdxyzQ(:,:,:) = dAdxyzQ(:,:,:) + tmp(:,:,:)

    end subroutine eemdAdxyzTimesQEwald

    subroutine dAdxyzQReal(nat, ats, lat, sigma, q, eta, dAdxyzQ)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), sigma(nat), q(nat)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: dAdxyzQ(nat, 3, nat)

        real(dp) :: d(3), r, dlat(3), inv2eta, invsqrt2eta, interf, gamma, f(3), d2, cutoff
        integer :: i, j, k, iat, jat, n(3)

        inv2eta = 1._dp / (2._dp * eta)
        invsqrt2eta = 1._dp / (sqrt(2._dp) * eta)


        cutoff = getOptimalCutoffReal(eta, ewaldSummationPrecision)
        call getNcells(lat, cutoff, n)

        dAdxyzQ(:,:,:) = 0._dp
        do iat=1,nat
            do i=-n(1),n(1)
                do j=-n(2), n(2)
                    do k=-n(3), n(3)
                        dlat(:) = i * lat(:,1) + j * lat(:,2) + k * lat(:,3)
                        do jat=1,nat
                            if (i/=0 .or. j/=0 .or. k/=0 .or. iat/=jat) then
                                d(:) = ats(:,iat) - ats(:,jat) + dlat(:)
                                d2 = sum(d**2)
                                if (d2 > cutoff**2) cycle
                                r = sqrt(d2)
                                interf = erfc(r * invsqrt2eta)
                                if (sigma(1) > 0._dp) then
                                    gamma = sqrt(sigma(iat)**2 + sigma(jat)**2)
                                    interf = interf - erfc(r / (sqrt(2._dp) * gamma))
                                end if

                                if (sigma(1) <= 0._dp) then
                                    f(:) = d(:) / r**3 &
                                            * (2._dp * r * invsqrt2eta * exp(-invsqrt2eta**2 * r**2) / sqrt(PI) &
                                                    + interf) / 2._dp
                                else
                                    f(:) = d(:) / r**3 &
                                      * (2._dp * r * invsqrt2eta * exp(-invsqrt2eta**2 * r**2) / sqrt(PI) &
                                      - 2._dp * r / (sqrt(2._dp) * gamma) * exp(-1._dp / (2._dp * gamma**2) * r**2) / sqrt(pi) &
                                      + interf) / 2._dp
                                end if
                                ! f = dA_{iat,jat}/dxyz
                                dAdxyzQ(iat, :, iat) = dAdxyzQ(iat, :, iat) - f(:) * Q(jat)
                                dAdxyzQ(jat, :, jat) = dAdxyzQ(jat, :, jat) + f(:) * Q(iat)
                                dAdxyzQ(jat, :, iat) = dAdxyzQ(jat, :, iat) - f(:) * Q(iat)
                                dAdxyzQ(iat, :, jat) = dAdxyzQ(iat, :, jat) + f(:) * Q(jat)
                            end if
                        end do
                    end do
                end do
            end do
        end do

    end subroutine dAdxyzQReal

    subroutine dAdxyzQrecip(nat, ats, lat, q, eta, dAdxyzQ)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat), lat(3,3), q(nat)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: dAdxyzQ(nat, 3, nat)
        real(dp) :: reclat(3,3), dlat(3), V, r, kri, krj, cutoff
        real(dp) :: factor
        integer :: i, j, k, iat, jat, n(3)
        real(dp) :: dkri(3), dkrj(3)
        real(dp) :: dAdxi(3)

        call recLattice(lat, reclat)
        V = det3D(lat)

        cutoff = getOptimalCutoffRecip(eta, ewaldSummationPrecision)
        call getNcells(reclat, cutoff, n)

        dAdxyzQ(:,:,:) = 0._dp

        do i=-n(1),n(1)
            do j=-n(2),n(2)
                do k=-n(3),n(3)
                    if (i/=0 .or. j/=0 .or. k/=0) then
                        dlat(:) = i * reclat(:,1) + j * reclat(:,2) + k * reclat(:,3)
                        r = sum(dlat**2) ! r = k**2
                        if (r > cutoff**2) cycle
                        factor = exp(-eta**2 * r / 2._dp) / r
                        do iat=1,nat
                            kri = sum(dlat(:) * ats(:,iat))
                            dkri(:) = dlat(:)
                            do jat=1,nat
                                krj = sum(dlat(:) * ats(:,jat))
                                dkrj(:) = dlat(:)
                                dAdxi(:) = factor * (-sin(kri) * dkri * cos(krj) + cos(kri) * dkri * sin(krj))

                                dAdxyzQ(iat, :, iat) = dAdxyzQ(iat, :, iat) + dAdxi(:) * Q(jat)
                                dAdxyzQ(jat, :, iat) = dAdxyzQ(jat, :, iat) + dAdxi(:) * Q(iat)

                                dAdxyzQ(jat, :, jat) = dAdxyzQ(jat, :, jat) - dAdxi(:) * Q(iat)
                                dAdxyzQ(iat, :, jat) = dAdxyzQ(iat, :, jat) - dAdxi(:) * Q(jat)

                            end do
                        end do
                    end if
                end do
            end do
        end do

        dAdxyzQ = dAdxyzQ / V * 4._dp * PI / 2._dp
    end subroutine dAdxyzQrecip

end module ewaldSummation
