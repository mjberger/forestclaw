      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, 2)

      do 30 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
            aux(i,j,1) = 0.0
   20       continue
   30       continue

       return
       end
