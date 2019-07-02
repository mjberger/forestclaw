      subroutine airburst_setprob(rho, bulk)
      implicit none
      double precision rho, bulk, cc, zz
      double precision rho_com, bulk_com, cc_com, zz_com, grav_com
      common /cparam/ rho_com,bulk_com,cc_com,zz_com, grav_com

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # These need to be assigned from user options
      rho_com = rho
      bulk_com = bulk
      cc = sqrt(bulk/rho)
      zz = rho*cc
      cc_com = cc
      zz_com = zz
      grav_com = 9.81d0

      return
      end
