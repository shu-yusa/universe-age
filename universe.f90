!======================================================================!
!     Title  : univ.f90                                                !
!     Author : Yusa Shusaku                                            !
!     Date   : 2008-12-16-Tue ~ 2008-12-17-Wed                         !
!     Last modified : 2008-12-18-Thu                                   !
!                                                                      !
!     A program which solves the evolution of the universe for several !
!     cosmological parameter sets.                                     !
!     The evolution of the universe is obtained by performing a        !
!     numerical integration with respect to scale factor 'a'.          !
!     In this program, we use Simpson's formula (with iteration) for   !
!     the integration. The results are accomodated to files and we     !
!     can plot them using, for example, gnuplot(this program produces  !
!     a file for gnuplot).                                             !
!                                                                      !
!     *** Structure of this program ***                                !
!                                                                      !
!     module      Com_var                                              !
!     module      Cosm_parameter                                       !
!     module      Cosm_parameter2                                      !
!     program     main                                                 !
!     function    Hub                                                  !
!     subroutine  Simpson                                              !
!     subroutine  Gnu                                                  !
!     subroutine  Title                                                !
!     subroutine  Show_Integral                                        !
!     subroutine  Show_Results                                         !
!                                                                      !
!======================================================================!
      module Com_var
!----------------------------------------------------------------------!
!     Definition of global constants.                                  !
!     kmax .... Maximum number of iteration of integration.            !
!     num_para .... The number of parameter sets.                      !
!     rad  .... 1 radian.                                              !
!----------------------------------------------------------------------!
      implicit none
      integer, parameter :: kmax=100, num_para=4
      real(8), parameter  :: epsr=1.0d-10
      real(8), parameter  :: rad=3.141592653589793d0/180.0d0
      end module
!======================================================================!
      module Cosm_parameter
!----------------------------------------------------------------------!
!     Definition of constants related to cosmology.                    !
!----------------------------------------------------------------------!
      implicit none
      real(8), parameter :: H0 = 1.0d0/13.6d0
      real(8), parameter :: a0 = 1.0d0
      end module
!======================================================================!
      module Cosm_parameter2
!----------------------------------------------------------------------!
!     Definition of derived data type.                                 !
!     This is used for cosmological parameter sets.                    !
!----------------------------------------------------------------------!
      implicit none
      type omega 
         real(8) :: Wm, Wr, Wk, Wl
      end type
      end module
!======================================================================!
      program main
!----------------------------------------------------------------------!
!     Main program.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : num_para
      use Cosm_parameter2 
      implicit none
      integer :: i
      real(8) :: t0(num_para)
      type(omega) :: W(num_para)
      character(len=20) :: Fname(num_para)

      call Title
      call Show_Integral

      W(1) = omega(1.0d0, 0.0d0,  0.0d0, 0.0d0)
      W(2) = omega(0.3d0, 0.0d0,  0.7d0, 0.0d0)
      W(3) = omega(0.3d0, 0.0d0,  0.0d0, 0.7d0)
      W(4) = omega(3.0d0, 0.0d0, -2.0d0, 0.0d0)

      do i=1, num_para
        Fname(i) = 'model'//char(48+i)//'.dat'
      end do
      
      call Solution(W, t0, Fname)
      call Gnu(7, Fname, 'gnu')
      call Show_Results(W, t0)

      end program
!======================================================================!
      subroutine  Solution(W, t0, Fname)
!----------------------------------------------------------------------!
!     This subroutine performs calculation for each parameter sets.    !
!     For parameter sets with non-positive curvature (Wk >= 0), we     !
!     perform numerical integration.                                   !
!     For parameter sets with positive curvature(Wk < 0) and without   !
!     dark energy(Wl=0), we use exact solution.                        !
!     Exact solution is also used for flat(Wk=0) universe without dark !
!     energy(Wl=0).                                                    !
!     We set the origin of the time to the value correspoding to a=1,  !
!     which is the present time.                                       !
!----------------------------------------------------------------------!
      use Com_var, only : num_para, rad
      use Cosm_parameter, only : H0
      use Cosm_parameter2
      implicit none
      integer :: i, j
      integer, parameter :: imax=200
      real(8),  parameter :: da=0.02d0
      real(8), intent(out) :: t0(num_para)
      real(8) :: eta, a, c, t
      type(omega), intent(in) :: W(num_para)
      character(*), intent(in) :: Fname(num_para)

      open(7, file=Trim(Fname(1)))
      a = 0.0d0
      do i=1, imax
         t0(1) = 2.0d0 / (3.0d0 * H0 * sqrt(W(1)%Wm)) 
         write(7,*) t0(1) * (a ** (1.5d0) - 1.0d0), a
         a = a + da
      end do
      close(7)

      do j=2, 3
         open(j+6, file=Trim(Fname(j)))
         call Simpson(0.001d0, 1.0d0, W(j), t0(j))
         a = 0.002d0
         do i=1, imax
            call Simpson(0.001d0, a, W(j), t)
            write(j+6,*) t - t0(j), a
            a = a + da
         end do
         close(j+6)
      end do

      eta = asin(sqrt(- W(4)%Wk / W(4)%Wm))
      c = W(4)%Wm / H0 / (- W(4)%Wk) ** (1.5d0)
      t0(4) =  c * (eta - sin(eta) * cos(eta))

      open(10, file=Fname(4))
      eta = 0.0d0
      do i=1, 180
         t = c * (eta - sin(eta) * cos(eta)) - t0(4)
         a = - W(4)%Wm / W(4)%Wk * sin(eta) ** 2
         write(10,*) t, a
         eta = eta + rad
      end do
      close(10)

      return
      end subroutine
!======================================================================!
      function  Hub(W, a)  result(H)
!----------------------------------------------------------------------!
!     Definition of an integrand.                                      !
!----------------------------------------------------------------------!
      use Cosm_parameter, only : H0, a0
      use Cosm_parameter2 
      implicit none
      real(8), intent(in) :: a
      real(8) :: H, c
      type(omega), intent(in) :: W

      c = a0 / a
      H = ((W%Wr * c + W%Wm) * c + W%Wk) * c * c + W%Wl
      H = 1.0d0 / (a * H0 * sqrt(H))

      return
      end function
!======================================================================!
      subroutine Simpson(a, b, W, S)
!----------------------------------------------------------------------!
!     A subroutine which perform numerical integration by using        !
!     Simpson's formula. In this subroutine, we increase the division  !
!     points until the result converges.                               !
!----------------------------------------------------------------------!
      use Com_var, only : epsr, kmax
      use Cosm_parameter2
      implicit none
      integer :: j, k, N
      real(8), intent(in)  :: a, b
      real(8), intent(out) :: S
      real(8), external :: Hub
      type(omega), intent(in) :: W
      real(8) :: s1, s2, s4, ds, SS, h

      h = b - a
      s1 = Hub(W,a) + Hub(W,b) 
      s2 = 0.0d0
      s4 = 0.0d0
      S = 0.5d0 * h * s1
      k = 0
      N = 1
      
      do 
        k = k + 1
        N = N * 2
        h = 0.5d0 * h
        s2 = s2 + s4
        s4 = 0.0d0       
        do j=1, N-1, 2
          s4 = s4 + Hub(W,a+dble(j)*h)
        end do
        SS = S                                      
        S = (h / 3.0d0) * (s1 + 4.0d0 * s4 + 2.0d0 * s2)
        ds = abs(S - SS)
        if ( ds < epsr * abs(S) ) then
           exit  
        end if
        if ( k > kmax ) then 
           write(6,*) 'The calculation did not converge.'
           exit
        end if  
      end do

      return
      end subroutine
!======================================================================!
      subroutine  Gnu(Fnum, DatFName, GnuName)
!----------------------------------------------------------------------!
!     This subroutine produces a file for gnuplot.                     !
!----------------------------------------------------------------------!
      use Com_var, only : num_para
      implicit none
      integer, intent(in) :: Fnum
      character(*), intent(in) :: DatFName(num_para)
      character(*), intent(in) :: GnuName

      open(Fnum,file=GnuName)
      write(Fnum,*) 'set term postscript eps enhanced color'
      write(Fnum,*) 'set output "universe.eps"'
      write(Fnum,*) 'unset key'
      write(Fnum,*) 'set size 0.6,0.6'
!      write(Fnum,*) 'set grid '
      write(Fnum,*) 'set xlabel "t [Gyr]"'
      write(Fnum,*) 'set ylabel "a / a_0"'
      write(Fnum,*) 'set xrange [-20:50]'
      write(Fnum,*) 'set yrange [0:4]'
      write(Fnum,*) 'set label "(K>0,{/Symbol L}=0)" at -2,0.5'
      write(Fnum,*) 'set label "(K=0,{/Symbol L}=0)" at 32,2.5'
      write(Fnum,*) 'set label "(K<0,{/Symbol L}=0)" at 22,3.6'
      write(Fnum,*) 'set label "(K=0,{/Symbol L}>0)" at  4,3.5'
      write(Fnum,*) 'set arrow from 25,3.4 to 27,3'
      write(Fnum,*) 'pl "'//Trim(DatFname(1))//'" w l, "'              &
     &                    //Trim(DatFname(2))//'" w l, "'              &
     &                    //Trim(DatFname(3))//'" w l, "'              &
     &                    //Trim(datFname(4))//'" w l'
      write(Fnum,*) 'set term x11'
      close(Fnum)

      return
      end subroutine
!======================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '            *****************************'
      write(6,*) '              Expansion of the Universe  '
      write(6,*) '            *****************************'
      write(6,*)

      return
      end subroutine
!======================================================================!
      subroutine Show_Integral
      implicit none

      write(6,*)'|****************************************************|'
      write(6,*)'|                                                    |'
      write(6,*)'| An equation which we solve(integrate) :            |'
      write(6,*)'|                                                    |'
      write(6,*)'|       a                                            |'
      write(6,*)'|      /                                             |'
      write(6,*)'|  1   [ da                   1                      |'
      write(6,*)'| ---  | --- ----------------------------------- = t |'
      write(6,*)'|  H0  ]  a  sqrt[Wm*c^3 + Wr*c^4 + Wk*c^2 + Wl]     |'
      write(6,*)'|      /                                             |'
      write(6,*)'|      0                                             |'
      write(6,*)'|                                                    |'
      write(6,*)'| where c = a0 / a .                                 |'
      write(6,*)'|****************************************************|'
      write(6,*) 
     
      return
      end subroutine
!======================================================================!
      subroutine  Show_Results(W, t0)
      use Com_var, only : num_para
      use Cosm_parameter2
      implicit none
      integer :: i
      real(8), intent(in) :: t0(num_para)
      type(omega), intent(in) :: W(num_para)
      character(len=50) :: Q
      parameter(Q='(2x,a,f4.1,a,f4.1,a,f4.1,a,5x,f6.2," Gyr",8x,a)')

      write(6,*) ' Fixed Cosmological Parameters :'
      write(6,*) ' H0 = 72 km/s/Mpc = (13.6Gyr)^-1 '
      write(6,*) ' Wr = 0'
      write(6,*)
      write(6,*) ' Results :'
      write(6,*) ' |********************************************|'
      write(6,*) ' |  (  Wm,  Wk,  Wl)  |  Age of the universe  |'
      write(6,*) ' |--------------------------------------------|'
      do i=1, num_para
       write(6,Q)'|  (',W(i)%Wm,',',W(i)%Wk,',',W(i)%Wl,')  |',t0(i),'|'
      end do
      write(6,*) ' |********************************************|'
      write(6,*)

      return
      end subroutine
!======================================================================!
