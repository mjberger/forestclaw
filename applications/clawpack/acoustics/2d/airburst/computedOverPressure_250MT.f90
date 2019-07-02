!
! ==========================================================================
! converted from MJA C program with parameters than model Chelyabinsk explosion:  520kt @ 29km alt
!
!        return overpressure (in percent of sea level STD ATM) as a function
!        of rad = radius from ground zero (km) and time in seconds
!        ...based on Friedlander wave, but scaled and tuned from actual runs.
!
double precision function  computedOverPressure(dist_in_km,time)

   implicit none

   ! Arguments
   real(kind=8), intent(in) :: dist_in_km, time

   ! local arguments - could be params, but for testing make them variables
   double precision ::  width, thick, speed
   double precision :: c, currentRadius, g
   double precision :: maxAmp1, maxAmp2, maxAmp3, maxAmp4

    !                    ...set consts for 250Mt TNT blast at 10km altitude
    maxAmp1  = 192.d0 !from MJA sims  ! *10 of Tunguska sized object 
    maxAmp2  = 192.d0 !from MJA sims  ! *10 of Tunguska sized object 
    maxAmp3  =  10.d0  !from MJA sims  ! *10 of Tunguska sized object 
    maxAmp4  = 68.d0  !from MJA sims  ! *10 of Tunguska sized object 

    width   = 30.d0
    thick   = 17.d0 
    speed   = 0.3915d0
    c       =  width/2.35482d0

    currentRadius    =  time * speed
   
    ! ...pulse -- functional fit from various blast simulations 
    ! mja model computes %overpressure, so need to convert to units
    ! so 6% becomes .06*101300 , done above in calling program
    if ( dist_in_km <= currentRadius) then
       g    =  0.d0 + maxAmp1 * exp(-4.0d0*(currentRadius/c        )**2)  ! /* ...Gaussian envelope */
       g    =  g    + maxAmp2 * exp(-2.0d0*(currentRadius/(1.8d0*c))**2)  ! /* ...Gaussian envelope */
       g    =  g    + maxAmp3 * exp(-2.0d0*(currentRadius/(5.d0*c))**8)   ! /* ...Gaussian envelope */
       g    =  g    + maxAmp4 * exp(-2.0d0*(currentRadius/(10.d0*c)))     ! /* ...Gaussian envelope */

       computedOverPressure = g*exp(-3.8d0*(currentRadius-dist_in_km)/(2.d0*thick)) *  &
                            (1.0d0 - 6.0d0*(currentRadius-dist_in_km)/(2.d0*thick))

    else
       computedOverPressure = 0.d0
    endif

end function computedOverPressure
