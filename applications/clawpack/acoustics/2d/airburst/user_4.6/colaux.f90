!
! --------------------------------------------------------
!
      subroutine colaux(aux,maux,mitot,mjtot,xlow,hx,time)

!     output the one row of aux (height) for plotting 
!     output x,height pairs at a given time,
!     this routine assumes it is always the top row not counting
!     ghost cells that needs to be printed. (So set max1d so that
!     there is no refinement in y

      !use amr_module, only :  hunit
      implicit double precision (a-h,o-z)


      !dimension aux(maux,mitot,mjtot)
      real (kind=8), intent(in) ::  aux(mitot,mjtot,maux)
      real (kind=8)::  auxout
      integer :: hunit
!      common /cparam/ rho,bulk,cc,zz,blastx_center,speed,      &
!                     ambient_pressure,gravity

      ambient_pressure = 101300.d0
      blastx_center = 0.d0
      speed = 391.5d0
                     
      open(555,file='fort.555',position='append');

      blastLoc = blastx_center + speed*time
      hunit = 555

      do i = 3, mitot-2
         x = xlow + (i-3+.5d0)*hx
         dist_in_km = (x - blastLoc) / 1000.d0   !convert m to km
!         externalPress = computedOverPressure(dist_in_km,time)/100.d0
!         ext = externalPress
         ext = aux(i,mjtot-3,1)  ! to make sure using real one
         auxout = aux(i,mjtot-2,1)  ! <-- last real col, not ghost cell
         if (abs(auxout) .lt. 1.d-20) auxout = 0.d0
         if (ext .lt. 1.d-20) ext = 0.d0
         write(hunit,100) time,x,auxout,ext
 !        write(hunit+1,100) x,ext
 100     format(4e15.7)
      end do

      close(555)
      
      return
      end subroutine colaux


