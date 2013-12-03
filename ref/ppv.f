      real function ppv(T)

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c     Calculate saturation water vapour pressure
c
c!Input Parameters:
c    T      Temperature (K)
c
c!Output Parameters:
c    ppv    Saturation vapour pressure (mb)
c
c!Revision History:
c 06/04 Collection 5  Marco Matricardi. ECMWF.  Original version.
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save 

      real     T                 ! Temperature [K] 
      real     esi,esw
      real     ti,t00,e00

c     ------------------------------------------------------------------
c     1. Constants
c     ------------------------------------------------------------------
      e00 = 611.21
      t00 = 273.16
      ti  = t00 - 23.

c     ------------------------------------------------------------------
c     2. Saturation vapour pressure
c     ------------------------------------------------------------------

      esw = e00*exp(17.502*(t-t00)/(t-32.19))
      esi = e00*exp(22.587*(t-t00)/(t+0.7  ))

      if(T.gt.t00)then
        ppv = esw                                          ! Water phase 
      elseif(t.gt.ti.and.t.le.t00)then
        ppv = esi+(esw-esi)*((t-ti)/(t00-ti))**2           ! Mixed phase
      elseif(t.le.ti)then
        ppv = esi                                          ! Ice phase 
      endif

      ppv=ppv/100.                     ! Conversion from [Pascal] to [mb]

      return
      end
