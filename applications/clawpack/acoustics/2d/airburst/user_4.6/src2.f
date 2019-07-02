      subroutine clawpack46_src2(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux,t,dt)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
c
c     # source terms for cylindrical symmetry in 2d Euler equations
c     # about y = 0, so radius = y
c     # y value at cell center is stored in aux(i,j,1)
c
      double precision rho_com, bulk_com, cc_com, zz_com, grav_com
      common /cparam/ rho_com,bulk_com,cc_com,zz_com, grav_com
c
c

      do 10 i=1,mx
       do 10 j=1,my

         q(i,j,3) = q(i,j,3) - dt * q(i,j,1) * grav_com / bulk_com
   10    continue

      return
      end
