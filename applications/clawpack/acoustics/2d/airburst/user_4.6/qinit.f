      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c
c     # Set initial conditions for q.
c     # Acoustics with smooth radially symmetric profile to test accuracy
c
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c
       pi = 4.d0*datan(1.d0)
       width = 500d0

       do 20 i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             ycell = ylower + (j-0.5d0)*dy
             q(i,j,1) = 0.d0
             q(i,j,2) = 0.d0
             q(i,j,3) = 0.d0
  20         continue
       return
       end
