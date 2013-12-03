      subroutine extend_profile( ni, pres, temp, mixr, lat, tf, wf )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Interpolate the NCEP GDAS1 temperature and moisture profiles to
c    the 101 standard pressure levels required for CIMSS transmittance
c    models.
c
c!Input Parameters:
c    NI      Number of input GDAS1 pressure levels
c    PRES    GDAS1 pressure levels (hPa) from maximum to minimum
c    TEMP    GDAS1 temperature (K) at PRES levels
c    MIXR    GDAS1 water vapor mixing ratio (g/kg) at PRES levels
c    LAT     Latitude (deg) of GDAS1 data
c
c!Output Parameters:
c    TF      Temperature (K) at 101 standard levels
c    WF      Water vapor mixing ratio (g/kg) at 101 standard levels
c            (NOTE: see variable pref for pressures at 101 standard levels)
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none

c ... arguments
            
      integer ni
      real pres( ni ), temp( ni ), mixr( ni ), lat, tf(101), wf(101)
      
c ... local variables

      integer k, nf
      parameter ( nf = 101 )
      real pref( nf ), lnpf( nf ), pi( 1000 ), lnpi( 1000 ), ti( 1000 ),
     &  wi( 1000 )

c ... data statements

c ... 101 standard CIMSS pressure levels (hPa)
      data pref    / 0.0050,    0.0161,    0.0384,    0.0769,    0.1370,
     +    0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,
     +    1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,
     +    5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,
     +   14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,
     +   29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,
     +   51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,
     +   83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,
     +  125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,
     +  180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,
     +  247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,
     +  328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,
     +  424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,
     +  535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,
     +  661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,
     +  802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,
     +  958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/


c ... Reverse the order of the input p, t, w arrays
      
      do k = 1, ni
        pi( k ) = pres( ni - k + 1 )
        ti( k ) = temp( ni - k + 1 )
        wi( k ) = mixr( ni - k + 1 )
      end do

c ... Convert pressure to ln(pressure)

      do k = 1, ni
        lnpi( k ) = log( pi( k ) )
      end do
      do k = 1, nf
        lnpf( k ) = log( pref( k ) )
      end do
      
c ... Interpolate to the 101 standard levels

      call interp( ni, lnpi, ti, nf, lnpf, tf )
      call interp( ni, lnpi, wi, nf, lnpf, wf )

c ... Extend temperatures from top of input profile to
c ... top of output profile
      
      do k = 1, nf
        if ( pref( k ) .lt. pi( 1 ) ) tf( k ) = -1.0
      end do
      call extem101( tf, lat )
      
      end
