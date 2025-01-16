!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Giovanni Bussi
! Adapted to QE by Andrea Ferretti & Layla Martin Samos
!
!----------------------------------
  MODULE coulomb_vcut_module
  !----------------------------------
  !
  IMPLICIT NONE
  PRIVATE

  !
  ! general purpose parameters
  !
  INTEGER,  PARAMETER :: DP=KIND(1.0d0)
  REAL(DP), PARAMETER :: PI     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: TPI    = 2.0_DP * pi
  REAL(DP), PARAMETER :: FPI    = 4.0_DP * pi
  REAL(DP), PARAMETER :: e2     = 2.0_DP
  REAL(DP), PARAMETER :: eps6   = 1.0E-6_DP
 
  !
  ! definitions
  ! 
  TYPE vcut_type
      REAL(DP)          :: a(3,3)
      REAL(DP)          :: b(3,3)
      REAL(DP)          :: a_omega
      REAL(DP)          :: b_omega
      REAL(DP), POINTER :: corrected(:,:,:)
      REAL(DP)          :: cutoff
      LOGICAL           :: orthorombic
  END TYPE vcut_type

  !
  PUBLIC :: vcut_type
  PUBLIC :: vcut_init
  PUBLIC :: vcut_get
  PUBLIC :: vcut_spheric_get
  PUBLIC :: vcut_destroy
  PUBLIC :: vcut_info

CONTAINS

!------------------------------------------
  SUBROUTINE vcut_init(vcut,a,cutoff)
  !------------------------------------------
  !
  TYPE(vcut_type),   INTENT(OUT) :: vcut
  REAL(DP),          INTENT(IN)  :: a(3,3)
  REAL(DP),          INTENT(IN)  :: cutoff

  INTEGER      :: n1,n2,n3
  INTEGER      :: i1,i2,i3
  INTEGER      :: ierr
  REAL(DP)     :: q(3)
  CHARACTER(9) :: subname='vcut_init'     
  REAL(DP)     :: mod2a(3)

  vcut%cutoff=cutoff

  vcut%a=a
  vcut%b= TPI * transpose(num_inverse(vcut%a))
  vcut%b_omega=num_determinant(vcut%b)
  vcut%a_omega=num_determinant(vcut%a)

  ! automatically finds whether the cell is orthorombic or not
  vcut%orthorombic=.false.
  !
  mod2a=sqrt(sum(vcut%a**2,1))
  if(abs(sum(vcut%a(:,1)*vcut%a(:,2)))/(mod2a(1)*mod2a(2))<eps6 .and. &
     abs(sum(vcut%a(:,2)*vcut%a(:,3)))/(mod2a(2)*mod2a(3))<eps6 .and. &
     abs(sum(vcut%a(:,3)*vcut%a(:,1)))/(mod2a(3)*mod2a(1))<eps6) vcut%orthorombic=.true.
  !
  if (.not.vcut%orthorombic) call errore(subname,"'vcut' Coulomb cutoff with non-orthogonal axis untested",1)

  n1=ceiling(vcut%cutoff*sqrt(sum(vcut%a(1,:)**2))/(2.0*pi))
  n2=ceiling(vcut%cutoff*sqrt(sum(vcut%a(2,:)**2))/(2.0*pi))
  n3=ceiling(vcut%cutoff*sqrt(sum(vcut%a(3,:)**2))/(2.0*pi))

  ALLOCATE(vcut%corrected(-n1:n1,-n2:n2,-n3:n3), STAT=ierr)
  IF ( ierr/=0 ) CALL errore(subname,'allocating cvut%corrected',ABS(ierr))
  !
  vcut%corrected=0.0

  !
  ! define the Fourier component of the modified Coulomb potential
  !
  DO i1=-n1,n1
    DO i2=-n2,n2
      DO i3=-n3,n3
        !
        q = MATMUL(vcut%b,(/i1,i2,i3/)) 
        !
        IF( SUM(q**2) > vcut%cutoff**2 ) CYCLE
        !
        vcut%corrected(i1,i2,i3) = &
             vcut_formula(q,vcut%a,vcut%b,vcut%a_omega,vcut%orthorombic)
        !
      ENDDO
    ENDDO
  ENDDO
  !
END SUBROUTINE vcut_init

!------------------------------------------
  SUBROUTINE vcut_info(iun, vcut)
  !------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN) :: iun
  TYPE(vcut_type), INTENT(IN) :: vcut
  !
  INTEGER :: i, n(3)
  !
  IF ( ASSOCIATED( vcut%corrected ) ) THEN
     !
     DO i = 1, 3
        n(i) = ( SIZE( vcut%corrected, i) -1 ) / 2
     ENDDO
     !
     WRITE(iun, "(  2x,'Cutoff: ',f6.2,4x,'  n grid: ',3i4,/)") vcut%cutoff, n(:)
     !
  ENDIF
  !
END SUBROUTINE vcut_info

!------------------------------------------
  SUBROUTINE vcut_destroy(vcut)
  !------------------------------------------
  !
  TYPE(vcut_type), INTENT(INOUT) :: vcut
  INTEGER :: ierr
  !
  DEALLOCATE(vcut%corrected, STAT=ierr)
  IF ( ierr/=0 ) CALL errore('vcut_destroy','deallocating vcut',ABS(ierr))
  !
END SUBROUTINE vcut_destroy

!------------------------------------------
  FUNCTION vcut_get(vcut,q) RESULT(res)
  !------------------------------------------
  !
  TYPE(vcut_type), INTENT(IN) :: vcut
  REAL(DP),        INTENT(IN) :: q(3)
  REAL(DP)                    :: res
  !
  REAL(DP)     :: i_real(3)
  INTEGER      :: i(3)
  CHARACTER(8) :: subname='vcut_get'
  !
  i_real=(MATMUL(TRANSPOSE(vcut%a),q))/ TPI
  i=NINT(i_real)
  !
  ! internal check
  IF( SUM( (i-i_real)**2 ) > eps6 ) &
     CALL errore(subname,'q vector out of the grid',10)
  !
  IF( SUM(q**2) > vcut%cutoff**2 ) THEN
     !
     ! usual form of Coulomb potential
     res = FPI * e2 / SUM(q**2) 
     !
  ELSE
     !
     IF( i(1)>ubound(vcut%corrected,1) .OR. i(1)<lbound(vcut%corrected,1) .OR. &
         i(2)>ubound(vcut%corrected,2) .OR. i(2)<lbound(vcut%corrected,2) .OR. &
         i(3)>ubound(vcut%corrected,3) .OR. i(3)<lbound(vcut%corrected,3)) THEN
         CALL errore(subname,'index out of bound', 10) 
     ENDIF
     !
     res=vcut%corrected(i(1),i(2),i(3))
     !
  ENDIF
  !
END FUNCTION vcut_get


!------------------------------------------
subroutine wofz(z, result_value)
  ! wofz computes the Faddeeva function w(z) for complex z,
  ! the funcion used here is derived from the plasma dispersion function z(x) = i*sqrt(pi)*w(z),
  ! the following plasma dispersion function is a multi-pole approximation by Huasheng Xie in 
  ! article "AIP Advances 1 July 2024; 14 (7): 075007 ", https://github.com/hsxie/gpdf/tree/main
  ! With a error of 1e-13 for J=24, which is enough for a double precision calculation here
  !------------------------------------------
  complex(dp), intent(in) :: z
  complex(dp) :: z_conj = cmplx(0.0_dp,0.0_dp, dp)
  complex(dp) :: b(12)
  complex(dp) :: c(12)
  integer    :: n
  complex(dp) :: result_value

  b( 1) = cmplx( -579.77656932346560644_dp , -844.01436313629880827_dp                  ,dp)
  b( 2) = cmplx( -179.52530851977905732_dp , -86.660002027244731382_dp                  ,dp)
  b( 3) = cmplx( -52.107235029274485215_dp , 453.3246806707749413_dp                    ,dp)
  b( 4) = cmplx( -2.1607927691932962178_dp , 0.63681255371973499384_dp                  ,dp)
  b( 5) = cmplx( -0.018283386874895507814_dp , -0.21941582055233427677_dp               ,dp)
  b( 6) = cmplx( -0.00006819511737162705016_dp , 0.00032026091897256872621_dp           ,dp)
  b( 7) = cmplx( -0.0000028986123310445793648_dp , -0.00000099510625011385493369_dp     ,dp)
  b( 8) = cmplx(  0.0000000023382228949223867744_dp , -0.0000000040404517369565098657_dp,dp)
  b( 9) = cmplx(  0.01221466589423530596_dp , 0.00097890737323377354166_dp              ,dp)
  b(10) = cmplx(  7.3718296773233126912_dp , -12.575687057120635407_dp                  ,dp)
  b(11) = cmplx( 44.078424019374375065_dp , -46.322124026599601416_dp                   ,dp)
  b(12) = cmplx( 761.62579175738689742_dp , 185.11797721443392707_dp                    ,dp)

  c( 1) = cmplx(   0.16167711630587375808393823760988_dp , -2.9424665391729649010502939606152_dp, dp)
  c( 2) = cmplx(   1.15091358764935672445993980434790_dp , -2.8745542965490153159866506667543_dp, dp)
  c( 3) = cmplx(   0.81513635269214329286824152984179_dp , -2.9085569383176322446978082849749_dp, dp)
  c( 4) = cmplx(   2.23629505890417241107360738208440_dp , -2.7033607074680388479084431872604_dp, dp)
  c( 5) = cmplx(   2.64035613134040415412304948466250_dp , -2.6228400297078984516779261304916_dp, dp)
  c( 6) = cmplx(   3.56204974511970566578349904839670_dp , -2.4245607245823420555878190731282_dp, dp)
  c( 7) = cmplx(   4.11692512571067539307286087375100_dp , -2.3036541720854573608940600179944_dp, dp)
  c( 8) = cmplx(   4.80341174933603179331098307177070_dp , -2.1592490859689535412501218722927_dp, dp)
  c( 9) = cmplx(   3.07789223492465673164827504614580_dp , -2.5301774598854448463007864644617_dp, dp)
  c(10) = cmplx( - 1.85720886352407650035610904791930_dp , -2.7720571884094886583775397071469_dp, dp)
  c(11) = cmplx(   1.49698813224668933803966639021490_dp , -2.8290855580900544693059801078858_dp, dp)
  c(12) = cmplx( - 0.48636891219330428093331493852099_dp , -2.9311741817223824196339069754696_dp, dp)
	
	

	z_conj = conjg(z)
	IF (aimag(z) >= 0._dp ) THEN
	  do n = 1, 12
        result_value = result_value + b(n)/(z-c(n)) + conjg(b(n))/(z+conjg(c(n)))
      end do
	ELSE
	  do n = 1, 12
        result_value = result_value + b(n)/(z_conj-c(n)) + conjg(b(n))/(z_conj+conjg(c(n)))
      end do
	  result_value = conjg(result_value)
	  result_value = result_value + cmplx(0.0_dp,2.0_dp, dp)*sqrt(pi)*exp(-z**2)
	ENDIF
    result_value = result_value / cmplx(0.0_dp,1.0_dp, dp) / sqrt(pi)
end subroutine wofz


subroutine erf_complex(z, result_value)
  ! erf_complex computes the error function for complex z
  ! please refer to http://ab-initio.mit.edu/faddeeva/
  implicit none
  INTEGER, PARAMETER :: dp = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: PI     = 3.14159265358979323846_dp
  complex(dp), intent(in)  :: z
  complex(dp)  :: tmp
  complex(dp)  :: result_value


    IF (real(z) >= 0._dp ) THEN
      tmp = cmplx(0.0_dp,1.0_dp, dp)*z
      call wofz(tmp, result_value)
	  result_value = cmplx(1.0_dp,0.0_dp, dp) - exp(-z**2) * result_value
	ELSE
      tmp = cmplx(0.0_dp,-1.0_dp, dp)*z
      call wofz(tmp, result_value)
	  result_value = exp(-z**2) * result_value-cmplx(1.0_dp,0.0_dp, dp)
	ENDIF
   
end subroutine erf_complex


!------------------------------------------
  FUNCTION vcut_spheric_get(vcut, q, exx_fraction, exx_fraction_lr, screen_parameter) RESULT(res)
  !------------------------------------------
  !
  TYPE(vcut_type), INTENT(IN) :: vcut
  REAL(DP),        INTENT(IN) :: q(3)
  REAL(DP)                    :: res 
  REAL(DP)                    :: res_sr
  REAL(DP),optional,INTENT(IN) :: exx_fraction
  REAL(DP),optional,INTENT(IN) :: exx_fraction_lr
  REAL(DP),optional,INTENT(IN) :: screen_parameter
  complex(dp) :: z, result_value

  !
  REAl(DP) :: a(3,3), Rcut, kg2 
  LOGICAL  :: limit
  !
  !

  a = vcut%a
  !
  Rcut=0.5*minval(sqrt(sum(a**2,1)))
  Rcut=Rcut-Rcut/50.0
  limit=.false.
  kg2=sum(q**2)
  result_value = cmplx(0.0_dp,0.0_dp, dp)
  if(kg2<eps6) then
    limit=.true.
  endif
  ! spherical cut method for bare Coulomb potential and its short-range part
  ! see e.g. J. Spencer and A. Alavi, Phys. Phys. Rev. B 77, 193110 (2008)
  ! please refer to https://www.vasp.at/wiki/index.php/Coulomb_singularity 
  ! for the details of the spherical cut method
  ! V_exx = exx_fraction*V_coulomb + exx_fraction_lr * V_coulomb_lr
  !       = (exx_fraction+exx_fraction_lr)*V_coulomb - exx_fraction_lr * V_coulomb_sr
  if(.not.limit) then
    ! full part
    res=TPI*e2/kg2*(1.0-cos(Rcut*sqrt(kg2)))
    if (exx_fraction_lr/=0._dp .and. screen_parameter > 0._dp) then
      z = cmplx(screen_parameter*Rcut, sqrt(kg2)/2.0_dp*screen_parameter, dp)
      call erf_complex(z, result_value)
      ! short-range part
      res_sr = TPI*e2/kg2*(1.0-cos(Rcut*sqrt(kg2))*erfc(screen_parameter*Rcut) - &
                           exp(-kg2/4.0_dp/screen_parameter**2)*REAL(result_value))
      res = (res*(exx_fraction+exx_fraction_lr) - res_sr*exx_fraction_lr)/exx_fraction
    end if
  else
    ! full part
    res=FPI*e2*Rcut**2/2.0
    if (exx_fraction_lr/=0._dp .and. screen_parameter > 0._dp) then
      ! short-range part
      res_sr = TPI*(Rcut**2*erfc(screen_parameter*Rcut) - &
                    Rcut*exp(-screen_parameter**2*Rcut**2)/(sqrt(pi)*screen_parameter) + &
                    erf(screen_parameter*Rcut)/2.0_dp/screen_parameter**2)
      res = (res*(exx_fraction+exx_fraction_lr) - res_sr*exx_fraction_lr)/exx_fraction
    end if
  endif
  !
END FUNCTION vcut_spheric_get

!---------------------------------------------------------
  FUNCTION vcut_formula(q,a,b,a_omega,orthorombic) result(res)
  !---------------------------------------------------------
  !
  ! Define the FT of the Coulomb potential according to the
  ! current lattice.
  !
  REAL(DP), INTENT(IN) :: q(3)
  REAL(DP), INTENT(IN) :: a(3,3)
  REAL(DP), INTENT(IN) :: b(3,3)
  REAL(DP), INTENT(IN) :: a_omega
  LOGICAL,  INTENT(IN) :: orthorombic
  REAL(DP)             :: res
  !
  real(dp) :: rwigner
  real(dp) :: sigma

  rwigner=0.5*sqrt(1.0/maxval(sum(b**2,1)))*2*pi

  !
  ! 3.0 is set to give a short range contribution inside the WS cell
  !
  sigma=3.0/rwigner

  ! compute longrange and shortrange contributions
  res=vcut_formula_longrange(q,a,b,a_omega,sigma,6.0D0,orthorombic) &
     +vcut_formula_shortrange(q,sigma)

END FUNCTION vcut_formula

!---------------------------------------------------------
  FUNCTION vcut_formula_longrange(q,a,b,a_omega,sigma,security,orthorombic) result(res)
  !---------------------------------------------------------
  ! compute the longrange contribution
  real(dp), intent(in) :: q(3)
  real(dp), intent(in) :: a(3,3)
  real(dp), intent(in) :: b(3,3)
  real(dp), intent(in) :: a_omega
  real(dp), intent(in) :: sigma
  real(dp), intent(in) :: security ! it determines the grid for the real-space sum; a reasonable value is 4.0
  logical,  intent(in) :: orthorombic
  real(dp)             :: res
  integer :: n1,n2,n3
  integer :: i1,i2,i3
  real(dp) :: d1,d2,d3,weight,factor
  real(dp) :: r(3),r2,modr
  logical :: n1_is_even,n1_is_odd
  real(dp) :: tmp, rtmp(3)
  logical, parameter :: shifted=.false.
  integer :: n1max
  real(dp) :: i1_real,i2_real,i3_real

  n1=security*sqrt(sum(a(:,1)**2))*sigma
  n2=security*sqrt(sum(a(:,2)**2))*sigma
  n3=security*sqrt(sum(a(:,3)**2))*sigma

  n1_is_even=(n1/2)*2==n1
  n1_is_odd=.not.n1_is_even

  d1=1.0/real(n1,dp)
  d2=1.0/real(n2,dp)
  d3=1.0/real(n3,dp)
  res=0.0
  weight=a_omega*d1*d2*d3
! the only symmetry which can be used for any value of q is inversion
! NON-SHIFTED:
!  if n1 is even: loop between 0 and n1/2, with weight=2.0 for all points except 0 and n1max
!  if n2 is odd:  loop between 0 and (n1+1)/2, with weight=2.0 for all points except 0
! SHIFTED:
!  if n1 is even: loop between 0 and n1/2-1, with weight=2.0 for all points
!  if n2 is odd:  loop between 0 and (n1+1)/2, with weight=2.0 for all points except n1max

  if(shifted)then
    if(n1_is_even) n1max=n1/2-1
    if(n1_is_odd)  n1max=(n1+1)/2-1
  else
    if(n1_is_even) n1max=n1/2
    if(n1_is_odd)  n1max=(n1+1)/2
  end if
  do i1=0,n1max
    factor=2.0
    if(shifted) then
      if(n1_is_odd .and. i1==n1max) factor=1.0
    else
      if(i1==0) factor=1.0
      if(n1_is_even .and. i1==n1max) factor=1.0
    end if
    i1_real=i1
    if(shifted) i1_real=i1_real+0.5
    do i2=0,n2-1
      i2_real=i2
      if(shifted) i2_real=i2_real+0.5
      do i3=0,n3-1
        i3_real=i3
        if(shifted) i3_real=i3_real+0.5
        rtmp=matmul(a,(/i1_real*d1,i2_real*d2,i3_real*d3/))
        r=vcut_minimal_image(a,b,rtmp,orthorombic)
        r2=sum(r**2)
        modr=sqrt(r2)
        if(modr*sigma<eps6) then
          tmp=e2*sqrt(2.0/pi)*sigma
        else
          tmp=e2*ERF(sigma*sqrt(0.5)*modr)/modr
        end if
        res=res+weight*factor*tmp*cos(sum(r*q))
      end do
    end do
  end do
 END FUNCTION vcut_formula_longrange

!---------------------------------------------------------
FUNCTION vcut_formula_shortrange(q,sigma) result(res)
  !---------------------------------------------------------
  real(dp), intent(in) :: q(3)
  real(dp), intent(in) :: sigma
  real(dp)             :: res
  if(sum(q**2/(sigma*sigma))<eps6) then
! analytic limit for small q
    res=e2*tpi/(sigma*sigma)
  else
    res=e2*fpi/sum(q**2)*(1-exp(-0.5*sum(q**2)/(sigma*sigma)))
  end if
END FUNCTION vcut_formula_shortrange

!---------------------------------------------------------
FUNCTION vcut_minimal_image(a,b,r,orthorombic) result(res)
  !---------------------------------------------------------
  real(dp), intent(in) :: a(3,3)
  real(dp), intent(in) :: b(3,3)
  real(dp), intent(in) :: r(3)
  logical,  intent(in) :: orthorombic
  real(dp)             :: res(3)
  real(dp) :: r_minimal(3)
  real(dp) :: r2_minimal
  real(dp) :: r_try(3)
  real(dp) :: r2_try
  real(dp) :: r_components(3)
  integer :: i1,i2,i3
  integer, parameter :: max_displacement=1
  if(orthorombic) then
! NINT ALGORITHM FOR ORTHOROMBIC CELL
    r_components=(matmul(transpose(b),r))/(2.0*pi)
    r_components=r_components-nint(r_components)
    r_minimal=matmul(a,r_components)
  else
! POOR MAN ALGORITHM FOR GENERIC CELL
    r_minimal=r
    r2_minimal=sum(r_minimal**2)
! loop over the possible neighbours
    do i1=-max_displacement,max_displacement
      do i2=-max_displacement,max_displacement
        do i3=-max_displacement,max_displacement
          if(i1==0 .and. i2==0 .and. i3==0) cycle
            r_try=r+matmul(a,(/i1,i2,i3/))
            r2_try=sum(r_try**2)
            if(r2_try<r2_minimal) then
              r2_minimal=r2_try
              r_minimal=r_try
            endif
        end do
      end do
    end do
  end if
  res=r_minimal
END FUNCTION vcut_minimal_image



!************************************
!** tools from sax

function num_inverse(a) result(inv)
  real(dp)              :: inv(0:2,0:2)
  real(dp), intent(in)  :: a(0:2,0:2)
  real(dp) :: tmp(0:2,0:2)
  real(dp) :: det
  real(dp),parameter :: eye3(3,3) = reshape((/ 1,0,0,0,1,0,0,0,1/),(/3,3/))
  integer i,j
  do i=0,2
    do j=0,2
      tmp(i,j) = a(modulo(i+1,3),modulo(j+1,3)) * a(modulo(i+2,3),modulo(j+2,3)) &
  &            - a(modulo(i+1,3),modulo(j+2,3)) * a(modulo(i+2,3),modulo(j+1,3))
    end do
  end do
  det = num_determinant(a)
  inv = transpose(tmp) / det
  if(sum((matmul(inv,a))**2-eye3) > 1d-5) then
    write(0,*) "AHIA",sum((matmul(inv,a)-eye3)**2)
    write(0,*) "A",a
    write(0,*) "inv",inv
    write(0,*)">>", matmul(inv,a)
    stop
  end if
end function num_inverse

function num_determinant(a) result(det)
  real(dp), intent(in) :: a(3,3)
  real(dp)             :: det
  det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) &
    - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3) - a(1,3)*a(2,2)*a(3,1)
end function num_determinant

!** end tools from sax
!************************************
END MODULE coulomb_vcut_module

