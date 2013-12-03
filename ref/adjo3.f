      subroutine adjo3(ozone, adj_ozone, mid_px, mid_ln)


      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'

      integer mid_px, mid_ln
      real ozone(plev), adj_ozone(plev)

      integer kk, jj, gdaspp(8), intpt, levc, m, n, strat_top,
     *        strat_bot, toplev, botlev
c     Define stratosphere boundaries (pressure levels).
      parameter ( strat_top = 21 )
      parameter ( strat_bot = 45 )

      real dp, do3,
     *     pp(plev),
     *     dptot,
     *     rlat, rlon
      
c ... Define 101 pressure levels.
      data pp      / 0.0050,    0.0161,    0.0384,    0.0769,    0.1370,
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

c     Define stratospheric interpolation points (pressure levels).
      data gdaspp /19, 21, 27, 30, 36, 40, 45, 47/

      if(debug .gt. 0) then
          write(h_output,'(''GDAS Total ozone: '',f10.4)') ozn1(mid_px,mid_ln)
          write(h_output,'(''mid_px, mid_ln: '',2i10)') mid_px, mid_ln
          write(h_output,'(''GDAS ozone profile: '',6f12.8)')
     &      o3prof1(6,mid_px,mid_ln),o3prof1(5,mid_px,mid_ln),o3prof1(4,mid_px,mid_ln),
     &      o3prof1(3,mid_px,mid_ln),o3prof1(2,mid_px,mid_ln),o3prof1(1,mid_px,mid_ln)

      end if

*****************************************************************************************

c     Copy original O3 profile to adjusted array.
      do kk = 1, plev
        adj_ozone(kk) = ozone(kk)
      enddo

*****************************************************************************************

c     Adjust climatology stratospheric ozone profile.
c     Replace climatological values with those from GDAS (6 values).
      adj_ozone(21) = o3prof1(6,mid_px,mid_ln)
      adj_ozone(27) = o3prof1(5,mid_px,mid_ln)
      adj_ozone(30) = o3prof1(4,mid_px,mid_ln)
      adj_ozone(36) = o3prof1(3,mid_px,mid_ln)
      adj_ozone(40) = o3prof1(2,mid_px,mid_ln)
      adj_ozone(45) = o3prof1(1,mid_px,mid_ln)

c     Interpolate stratospheric values.
      levc = 2
      intpt = 1
      toplev = strat_top - 2
      botlev = strat_bot + 2
      do kk = toplev, botlev

c       write(h_output,'(3i5,2f10.4)') kk,levc,intpt,ozone(gdaspp(levc)),
c    *          ozone(gdaspp(intpt))

        if(gdaspp(levc) .eq. kk) then
          do3 = adj_ozone(gdaspp(levc)) - adj_ozone(gdaspp(intpt)) 
          dptot = log(pp(gdaspp(levc))) - log(pp(gdaspp(intpt)))
c         write(h_output,'(3i5,4f10.4)') kk,levc,intpt,adj_ozone(gdaspp(levc)),
c    *          adj_ozone(gdaspp(intpt)),do3,dptot

          m = gdaspp(intpt) + 1
          n = gdaspp(levc) - 1
          do jj = m, n
            dp = log(pp(jj)) - log(pp(gdaspp(intpt)))
            adj_ozone(jj) = adj_ozone(gdaspp(intpt)) + (do3 * (dp / dptot))
c           write(h_output,'(''interp val: '',i5,4f10.4)') jj,ozone(jj),adj_ozone(jj),
c    *         dp,dp/dptot
          enddo

          levc = levc + 1
          intpt = intpt + 1

        end if

      enddo

*****************************************************************************************

      return
      end
