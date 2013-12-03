      subroutine get_co2cld( line, pixel, obs_clr, tcold, twarm, 
     &   met_date, pct, eca, qual_flag, conf_flag, co2_flag, bias_flag,
     &   hc_flag, ci_flag, os_top_flag, cldhgt_cat, nearnad_flag,
     &   cldhgtmet_flag, ncloudy, ngood_co2 )

c!F77------------------------------------------------------------------
c!Description: Compute cloud top pressure and effective cloud amount 
c              from MODIS IR data.
c
c!Input parameters:
c line                       Line number within a swath (1-10)
c pixel                      Pixel number within a 1km scan (1-1500)
c obs_clr                    Observed clear sky pixels
c                              0=none present
c                              1=present
c                              use only with mod06ct_proc_opt = 2
c tcold                      Obs cloud or cold brightness temperatures
c                              (K) for MODIS IR bands 36, 35, 34, 33, 31
c twarm                      Obs clear or warm brightness temperatures
c                              (K) for MODIS IR bands 36, 35, 34, 33, 31
c met_date                   Year, month, day, hour of gridded met data
c ncloudy                    Number of cloudy pixels
c ngood_co2                  Number of valid pixels
c
c The following in COMMON /MOD06_DATA/ are used:
c wprof1                     Water vapor mixing ratio profile (g/kg)
c tprof1                     Temperature profile (K)
c tem1                       Surface temperature (K)
c view                       Satellite zenith angle (deg)
c cs_bias_corr               clear-sky radiance bias correction 
c bad_value                  missing or bad data flag
c
c!Output parameters:
c pct                        Cloud top pressure (hPa)
c eca                        Cloud effective emissivity
c qual_flag                  Product quality flag (0=not useful, 1=useful)
c conf_flag                  Product confidence flag (0-3, 0=worst, 3=best)
c co2_flag                   Indicates whether or not co2 slicing was used
c                                   (0=no, 1 =yes)
c bias_flag                  Indicates whether or not zonal mean bias and instrument
c                              bias was used (0=no, 1=yes)
c hc_flag                    Indicates whether or not high cloud was
c                              detected (0=missing, 1=no, 2=yes)
c ci_flag                    Indicates whether or not cirrus cloud was
c                              detected (0=missing, 1=no, 2=yes)
c os_top_flag                Indicates over-shooting thunderstorm top (50S-50N only)
c cldhgt_cat                 Indicates lo, mid, hi clouds for QA output
c nearnad_flag               Indicates near-nadir view, VZA <= 32 deg. (1=yes, 2=no)
c cldhgtmet_flag             Indicates CTP retrieval method for QA output
c
c The following are output to COMMON /MOD06_DATA/:
c height_tropopause          Pressure of tropopause from temperature profile
c cloud_h_method             Number of the selected ratio solution
c cloudtop_pre_ir            Window cloud top height (mb)
c spec_cloud_forcing         Cloud-clear or cold-warm radiance for MODIS 
c                              IR bands 36, 35, 34, 33, 31 (mW/m2/sr/cm-1)
c cloudtop_pres              Final cloud top pressure (mb)
c cloudtop_height            Final cloud top height (meters)
c cloudtop_temp              Cloud top temperature from temperature profile
c cloudtop_eff_emi           Effective cloud amount (no units)
c cloudtop_pres_from_ratios  Cloud top pressure for each of 5 ratios
c
c!Revision history:
c  R. Frey 02-01-11
c  Changes for Collection 6
c  Implemented "top-down" method of final CTP choice for Aqua; solutons are
c  chosen in order of (36/35, 35/34, 34/33).
c  Restrict range of CTP retrievals appropriate to channel pair (36/35 <
c  450 hPa, 35/34 < 550 hPa, 34/33 < 650 hPa, 35/33 < 650 hPa).
c  Avoid CO2 slicing solutions in water clouds and IRW solutions in ice or
c  mixed phase clouds by use of "beta ratio" threshold (water surfaces only).
c  Use GDAS ozone profile data in stratosphere; merge with climatological
c  profiles currently in use.
c  Reduce NEDR thresholds for band selection in CO2-slicing algorithm.
c  Implement "spectral shift" in forward model calculations involving
c  involving bands 34-36 (Aqua only). Added routines modis_planck_shift and
c  modis_bright_shift.
c  Attempt to identify stratospheric clouds ("overshooting tops") by use of
c  13.6-11 um BTDs. Added flag 'os_top_flag'.
c  Implement alternate CTP algorithm for low clouds over oceans. Use clear
c  minus cloudy 11 um BT and precomputed lapse rates to generate cloud top
c  height, then convert to CTP. Added routine marine_locld_retrieval.
c  Output cloud (geopotential) heights corresponding to cloud top pressures.
c  Output subset of CTP, CTT, CTH, ECE, CF (cloud fraction) for "near-nadir"
c  views (VZA < 32 deg.).
c  Output new SDSs for L3 aggregation (solar zenith day and night, solar
c  azimuth day and night, sensor zenith day and night, sensor azimuth day
c  and night).
c  Added new QA flags for near-nadir, "over-shooting top", and cloud height
c  category (low, mid, and high clouds defined by 440mb and 680mb boundaries).
c  Removed checks for previously accessed temperature and moisture profiles. Code
c  now performs RT calculations on all profiles within the area covered by a
c  given granule.

c  R. Frey, H. Zhang 05-11-06
c  Changes for Terra only this time
c  Terra 8-day clear-sky zonal mean radiance for bands 31-36 and instrument bias
c  for bands 34-36 were added, adjusted rmin thresholds, no ratios 35/34 and
c  34/33 used, "Top-Down" method used for final CTP
c
c  R. Frey, H. Zhang 10-12-05
c  Aqua 8-day clear-sky zonal mean radiance for bands 31-36 and measurement/FM 
c  correction bias for bands 34-36 were added, zonal mean bias are
c  seperated for ocean, land/day, and land/night, for latitudes -60 to -90, biases 
c  were created from both day and night.
c  add check in tropopause pressure level determination; check for inflection point
c  in temperature profile.
c  co2_flag and bias_flag were added for tracking the frequencies of co2 slicing 
c  method and bias correction
c
c  R. Frey, H. Zhang 12-20-04
c  Use average cloud radiance from cloud pixels, cloud fraction used in IR
c  as ECA, also effect ECA from CO2 algorithm, no ratio 33/31 used,
c  also no ratio 35/33 for aqua. Took off the NEDR check at
c  11um for starting the co2 slicing algorithm, no low clouds quality control. 
c  selected bands used in finding the final co2 heights, band 33, 35 and 36 for
c  Terra, band 33 through 36 for Aqua.
c
c R. Frey  06-11-04
c Set Aqua "rmins" to 1.0, implemented 7.4 forward model, use more
c atmospheric levels in GDAS profiles as input, added surface emissivities,
c reduced number of forward model calculations
c 05-05-04 R. Frey, P. Hubanks, R. Hucek
c Implemented fix to Cirrus Flag
c R. Frey  10-17-03
c Implemented 101-level processing, climatological ozone profiles,
c atmospherically corrected window channel cloud pressures, check on
c cloud pressure functions, modified Aqua "rmins" for Aqua Collection 4
c delivery.  
c
c calls:
c function fm_modrad_emis     Compute MODIS IR TOA radiance given transmittance
c function modis_planck       Compute Planck radiance for MODIS IR bands
c function modis_planck_shift Compute Planck radiance for MODIS IR bands
c                             uses "shifted" central wave numbers for bands 27, 28, 34-36
c function modis_bright_shift Compute brightness temperature for MODIS IR bands
c                             uses "shifted" central wave numbers for bands 27, 28, 34-36
c subroutine tran_modisd101   Compute transmittance profile for MODIS IR bands
c subroutine clozo101         Computes ozone concentration profile from climatology
c subroutine height           Computes geopotential height
c subroutine adjo3            Adjusts stratospheric ozone profile using GDAS values
c subroutine assign_eco_emis  Assigns surface emissivities by frequency
c subroutine getiremis        Computes surface emissivity
c subroutine marine_locld_retrieval   Computes marine low cloud heights
c
c!Team-unique Header:
c
c!End------------------------------------------------------------------

      implicit none
      save

c----------------------------------------------------------------------

      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'
      include 'platform_name.inc'

c----------------------------------------------------------------------

      integer nb_wavelen
      parameter (nb_wavelen = 7)

c----------------------------------------------------------------------

c     Scalar arguments.
      real pct, eca
      integer qual_flag, conf_flag, ci_flag, hc_flag, line, pixel, obs_clr,
     &        ncloudy, ngood_co2, bias_flag, num_biasflag, co2_flag, os_top_flag,
     &        cldhgt_cat, nearnad_flag, cldhgtmet_flag

c----------------------------------------------------------------------

c     Array arguments.
      real tcold(ntbct), twarm(ntbct)
      integer met_date(4)

c----------------------------------------------------------------------

c     Local scalars
      real bot, top, db, fm, fm1, fm2, fmsav, sum, psfc, pmsl, tmin,
     &     tw, pfco2, ecawin, ecaco2, ptrp, pwin, tct, ts, theta, tss,
     &     rlon, rlat, emis_out, rho_out, emiswc, ppsfc, emisrw, rtmp1,
     &     rtmp2, zs, robs11, robs12, ne11, ne12, ne86, beta, emis12, 
     &     emis86, beta2, robs86, aqua_beta, terra_beta, beta_threshold,
     &     lapse_rate, pmlc
      integer id, is1, iw1, k, k1, k2, ll, ltrp, lwin, ip, mid_px,
     &        isp, imslp, lco2, ipco2, ngch,ipw, iout, jout, mid_ln,
     &        ii, jj, iisp, kban, jdet,
     &        met_month, lsf, landsea, init,
     &        llwin, nl, beta_lev, met_year, met_day, jday, lmin
      logical ok, neg, start
   
c----------------------------------------------------------------------

c     Local arrays.
      real amo(nsct), ra(plev,nbct), rclr(nbct), robs(nbct),
     &     emis(nb_wavelen), rho(nb_wavelen), freqemis(nb_wavelen),
     &     ozone(plev), pp(plev), tp(plev), ttpp(plev),
     &     taup(plev,ntbct), rmin(nbct), rwarm(nbct), rwcld(nsct),
     &     em_adj(nsct), delr(nbct), w(plev),
     &     bias(nbct), rmin_Terra(nbct), rmin_Aqua(nbct), rclr_s(plev),
     &     freq_terra(ntbct), freq_aqua(ntbct), sfc_emis(ntbct),
     &     aqua_rad_bias(nbct), terra_rad_bias(nbct), rad_bias(nbct),
     &     adj_ozone(plev), z(plev), rc1(nbct), rc2(nbct), rclr_s2(plev),
     &     ratio(nsct), freq(ntbct), rclr_s3(plev)
      integer lev(nsct), mbnds(ntbct), std_doy(12), leap_doy(12),
     &        krto(nsct), mcsbnds(nbct), kch(nsct,2)

c----------------------------------------------------------------------

c     External functions.
      real modis_bright_shift, modis_planck, modis_planck_shift, fm_modrad_emis
      external modis_bright_shift, modis_planck, modis_planck_shift, fm_modrad_emis

c----------------------------------------------------------------------

c     External subroutines.
      external adjo3, height, tran_modisd101, clozo101, getiremis,
     &         assign_eco_emis
c----------------------------------------------------------------------

c     Data statements.

c----------------------------------------------------------------------

c     Define 101 pressure levels.
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

c     Define indices of channel combinations used in the co2-slicing
c     method.
      data kch /1,2,2,3,2,3,4,4/

c     Initialize ozone profile.
      data ozone/plev*0.0/

c     Order of MODIS IR band numbers used
      data mbnds/36,35,34,33,31,29,32/      

c     Define wavenumbers of CO2 bands used (36, 35, 34, 33, 31).
      data freq_terra /7.055309E+02, 7.196677E+02, 7.317089E+02,
     &                 7.483224E+02, 9.081998E+02, 1.173198E+03,
     &                 8.315149E+02/
      data freq_aqua  /7.045020E+02, 7.190090E+02, 7.315760E+02,
     &                 7.482977E+02, 9.076808E+02, 1.169637E+03,
     &                 8.308397E+02/

c     Order of clear-sky radiance bias values used.
      data mcsbnds/5,4,3,2,1/

c     Ice cloud emissivity adjustments for each channel combination.
      data em_adj / 1.000, 1.000, 1.000, 1.000 /

c     Define noise level (in radiance units -- mw/m2/ster/cm-1) for each CO2 channel.
      data rmin_Aqua /-1.25, -1.0, -4.0, -4.0, -0.5/  
      data rmin_Terra /-1.0, -1.0, -100.0, -1.0, -0.5/  

c     Dectector number for transmittance calculations (0=averaged).
      data jdet /0/

c     Initialize 'init' flag
      data init /0/

c     Black cloud emissivity
      data emisrw /1.0/

c     Aqua radiometric biases (bands 31, 33-36 in that order).
      data aqua_rad_bias / 0.00, 0.00, 0.00, 0.00, 0.00 /
c     data terra_rad_bias / 0.00, 0.00, 0.66, 1.00, 1.50 /
      data terra_rad_bias / 0.00, 0.00, 0.00, 0.00, 0.00 /

c     Beta ratio thresholds
      data aqua_beta /0.95/
      data terra_beta /0.95/

      data std_doy  /0,31,59,90,120,151,181,212,243,273,304,334/
      data leap_doy /0,31,60,91,121,152,182,213,244,274,305,335/

c----------------------------------------------------------------------

c     Perform one time initializations (first call only).

      if(init .eq. 0) then

c       Check the platform name (it should be 'Terra' or 'Aqua')
        if (platform_name(1:5) .eq. 'Terra') then
          beta_threshold = terra_beta
          do k = 1, nbct
            rmin(k) = rmin_Terra(k)
            rad_bias(k) = terra_rad_bias(k)
          enddo 
          do k = 1, ntbct
            freq(k) = freq_terra(k)
          enddo
        else if (platform_name(1:4) .eq. 'Aqua') then
          beta_threshold = aqua_beta
          do k = 1, nbct
            rmin(k) = rmin_Aqua(k)
            rad_bias(k) = aqua_rad_bias(k)
          enddo 
          do k = 1, ntbct
            freq(k) = freq_aqua(k)
          enddo
        end if 

        init = 1

      end if

c----------------------------------------------------------------------

c     Initializations performed on every call.

      ngch = 0
      cldhgt_cat = 2
      cldhgtmet_flag = 5
      do id = 1,nsct
        krto(id) = 0
        lev(id) = 0
      enddo
      do k = 1, plev
        z(k) = 0.0
      enddo

c----------------------------------------------------------------------

c     Definitions for current 5x5.

c     Define output array indices for the current line and pixel.
      iout = (pixel / isamp) + 1
      jout = (line / isamp) + 1

c     Define mid-point of box in terms of line and pixel.
      mid_px = pixel + isamp / 2
      mid_ln = line + isamp / 2

c     Get satellite zenith angle.
      theta = view(mid_px,mid_ln)
      if(iout .gt. 675 .and. iout .lt. 679 .and. theta .lt. 0.0) then
        theta = 0.1
      end if

c     Get current latitude.
      rlat = lat1(mid_px,mid_ln)
      rlon = lon1(mid_px,mid_ln)

c     Get month, year, day-of-month from gridded met data.
      met_month = met_date(2)
      met_year = met_date(1)
      met_day = met_date(3)

c----------------------------------------------------------------------

c     Get climatological ozone profile.
      call clozo101(rlat, met_month, ozone)

c----------------------------------------------------------------------

c     Get surface pressure.
      psfc = pre1(mid_px,mid_ln)
c     Find level of surface pressure.
      do ll = 1,plev
        if(psfc .le. pp(ll)) then
          isp = ll
          go to 100
        end if
      enddo
 100  continue

c----------------------------------------------------------------------

c     Get temperature and water vapor mixing ratio profiles.
      do ii = 1,plev
        tp(ii) = tprof1(ii,mid_px,mid_ln)
        w(ii) = wprof1(ii,mid_px,mid_ln)
      end do

c----------------------------------------------------------------------

c     Get MSL pressure.
      pmsl = premsl(mid_px,mid_ln)

c     Find level of mean sea level pressure.
      do ll = 1,plev
        if(pmsl .le. pp(ll)) then
          imslp = ll
          go to 200
        end if
      enddo
 200  continue

c----------------------------------------------------------------------

c     Get geopotential height profile (km).
      zs = 0.0
      nl = imslp
      call height(pp,tp,w,zs,nl,z)

c----------------------------------------------------------------------

c     Get surface temperature.
      ts = tem1(mid_px,mid_ln)

c----------------------------------------------------------------------

c     Get land/sea tag.
      lsf = land1(mid_px,mid_ln)
      if(lsf .eq. 1 .or. lsf .eq. 2 .or. lsf .eq. 4) then
        landsea = 1
      else
        landsea = 0
      end if

c----------------------------------------------------------------------

c     Get surface emissivity.
      call assign_eco_emis(landsea,emis,rho,freqemis)
      do k = 1, ntbct
        call getiremis(nb_wavelen,freqemis,freq(k),emis,rho,emis_out,
     &                 rho_out)
        sfc_emis(k) = emis_out
      enddo

c----------------------------------------------------------------------

c     Adjust to ozone profile according to GDAS data.
      call adjo3(ozone, adj_ozone, mid_px, mid_ln)
 
c----------------------------------------------------------------------
 
c     Calculate day-of-year.
      if( mod(met_year,4) .eq. 0 ) then
        jday = leap_doy(met_month) + met_day
      else
        jday = std_doy(met_month) + met_day
      end if

c----------------------------------------------------------------------
 
c     Calculate transmittance profiles (101 level fast model).

c----------------------------------------------------------------------

      do k = 1, ntbct
        kban = mbnds(k)
        call tran_modisd101(met_year, jday, tp, w, adj_ozone, theta,
     &                      kban, jdet, taup(1,k))
      enddo

c----------------------------------------------------------------------
 
c     Write debug information for transmittance calculations.
c     if(debug .gt. 0) then
c       if(line .eq. 6 .and. pixel .eq. 676) then
c       if(line .eq. 6 .and. pixel .eq. 26) then
c       if(line .eq. 6 .and. pixel .eq. 1256) then
          write(*,'(1x)')
          write(*,'(2f10.3)') rlat, rlon
          write(*,'(''Profile Date: '',5i10)') (met_date(jj),jj=1,4), jday

          write(*,'(''GDAS ozone profile: '',6f12.8)')
     &      o3prof1(6,mid_px,mid_ln),o3prof1(5,mid_px,mid_ln),o3prof1(4,mid_px,mid_ln),
     &      o3prof1(3,mid_px,mid_ln),o3prof1(2,mid_px,mid_ln),o3prof1(1,mid_px,mid_ln)

          write(*,'(''Atm level, pressure, temperature, '',
     &      ''water vapor, ozone, transmittances (36,35,34,33,31,29,32)'')')
          do k=1,plev
            write(*,'(i4,12f10.4)') k,pp(k),tp(k),w(k),ozone(k), adj_ozone(k),
     &            (taup(k,jj),jj=1,ntbct)
          enddo
          read(*,*)
c       end if
c     end if

c----------------------------------------------------------------------

c     Find tropopause level.

c----------------------------------------------------------------------

c     Find level of coldest temperature between 100 mb and the surface.
      tmin = 350.0
      do ll = 44,isp
        if(tp(ll) .le. tmin) then
          tmin = tp(ll)
          lmin = ll
        end if
      enddo

c     Look for inflection point in temperature profile.
c     If found do not change tropopause level. 
c     In isothermal conditions, choose the bottom level.
c     If temperature increases monotonically below 100 mb, do not   
c     adjust tropopause level.

      ltrp = -99
      do ll = lmin, 71
        if(tp(ll-1) .ge. tp(ll) .and. tp(ll+1) .gt. tp(ll)) then
          ltrp = ll
          go to 300
        end if
      enddo

 300  continue
      if(ltrp .eq. -99) ltrp = lmin
      ptrp = pp(ltrp)

c----------------------------------------------------------------------

c     Get IRW (11 micron) radiance profile, then convert to BT.
      Do  k = ltrp,isp
        if(k .eq. isp) then
          emiswc = sfc_emis(kwc)
          emis12 = sfc_emis(7)
          emis86 = sfc_emis(6)
          ppsfc = psfc
          tss = ts
        else
          emiswc = 1.0
          emis12 = 1.0
          emis86 = 1.0
          ppsfc = pp(k)
          tss = tp(k)
        end if
c       tss = tp(k)
        iisp = k
        rclr_s(k) = fm_modrad_emis(tss,emiswc,ppsfc,pp,tp,taup(1,kwc),
     &                           mbnds(kwc),iisp)
        TTPP(k) = modis_bright_shift(rclr_s(k),mbnds(kwc),0)

c       Also get 12 and 8.6 micron radiance profile.
        rclr_s2(k) = fm_modrad_emis(tss,emis12,ppsfc,pp,tp,taup(1,7),
     &                             mbnds(7),iisp)
        rclr_s3(k) = fm_modrad_emis(tss,emis86,ppsfc,pp,tp,taup(1,6),
     &                             mbnds(6),iisp)
      enddo
  
c----------------------------------------------------------------------

c     Compute the "IR window" cloud height.  Compare 11 micron BT to
c     atmospherically-corrected temperatures from radiance profile.

c----------------------------------------------------------------------

      tw = tcold(kwc)
  
      lwin = ltrp
      do ll = ltrp,isp
        if(ttpp(ll) .lt. tw) then
          lwin = ll
        end if
      enddo

c     Find which level tw is closest to and define beginning level at
c     which to begin comparison of LHS and RHS of co2-slicing equation.
      if(lwin .eq. isp) then
        iw1 = lwin - 1
      else if(lwin .eq. ltrp) then
        iw1 = lwin
      else if( abs(tw - ttpp(lwin)) .gt. abs(tw - ttpp(lwin + 1)) ) then
        lwin = lwin + 1
        iw1 = lwin - 1
      else
        iw1 = lwin
      end if

      lco2 = lwin

c----------------------------------------------------------------------

c     Define all output values for window cloud height retrieval.
      pwin = pp(lwin)
      ecawin = float(ncloudy) / ngood_co2 
      ipw = 6
c     Round cloud pressure to nearest 5 mb.
      if(pwin .ge. psfc) pwin = psfc
      pwin = (nint(pwin / 5.0)) * 5.0
      pct = pwin
      eca = ecawin
      ip = ipw
      tct = tp(lwin)

c     Cloud Height Method QA Flag
      cldhgtmet_flag = 6

c     Initial values of confidence and quality flags. May be reset later.
      qual_flag = 1
      conf_flag = 3

c----------------------------------------------------------------------

c     Fill all output arrays which pertain to the "window" cloud height
c     retrieval. Values may be overwritten below.
      cloud_h_method(iout,jout) = ip * 1.0
      cloudtop_pres(iout,jout) = pct
      cloudtop_height(iout,jout) = (nint((z(lwin)*1000.0) / 50.0)) * 50.0
      if(cloudtop_height(iout,jout) .lt. 0.0) cloudtop_height(iout,jout) = 0.0
      cloudtop_temp(iout,jout) = tct
      height_tropopause(iout,jout) = (nint(ptrp / 5.0)) * 5.0
      cloudtop_eff_emi(iout,jout) = eca
      cloudtop_pre_ir(iout,jout) = pwin

c----------------------------------------------------------------------

c     Compute Beta ratio at average level between IR window cloud pressure
c     and tropopause.
      beta_lev = (lco2 + ltrp) / 2
      robs11 = modis_planck(tcold(5),31,0)
      robs12 = modis_planck(tcold(7),32,0)
      robs86 = modis_planck(tcold(6),29,0)
      ne11 = (robs11 - rclr_s(isp)) / (rclr_s(beta_lev) - rclr_s(isp))
      ne12 = (robs12 - rclr_s2(isp)) / (rclr_s2(beta_lev) - rclr_s2(isp))
      ne86 = (robs86 - rclr_s3(isp)) / (rclr_s3(beta_lev) - rclr_s3(isp))
c     12/11 micron beta-ratio.
      if(ne11 .le. 1.0 .and. ne12 .le. 1.0) then
        beta = log(1.0 - ne12) / log(1.0 - ne11)
      else
        beta = 999.9
      end if
c     8.6/11 micron beta-ratio.
      if(ne11 .le. 1.0 .and. ne86 .le. 1.0) then
        beta2 = log(1.0 - ne86) / log(1.0 - ne11)
      else
        beta2 = 999.9
      end if

c----------------------------------------------------------------------

c     Write debug information about input data, "IR window" cloud
c     height retrieval, and beta ratio.
      if(debug .gt. 0) then
c      if(line .eq. 6 .and. pixel .eq. 676) then
c      if(line .eq. 6 .and. pixel .eq. 26) then
c      if(line .eq. 6 .and. pixel .eq. 1256) then
        write(h_output,'(1x)')
        write(h_output,'(''get_co2cld debug'')')
        write(h_output,'(''Line and pixel '',2i5)') line,pixel
        write(h_output,'(''Latitude, longitude, VZA, LST '',3f10.3,i5)') rlat,
     &      rlon, theta, landsea
        write(h_output,'(''Output array indices are '',2i10)')
     &      iout,jout
        write(h_output,'(''Mid-point of box is at '',2i10)')
     &      mid_ln,mid_px
        write(h_output,'(1x)')
        write(h_output,'(''Input data / window height retrieval'')')
        write(h_output,'(''Trop level and pressure, beta ratios '',i5,3f10.3)')
     &      ltrp, ptrp, beta, beta2
        write(h_output,'(''Surface level, temp, and pressure '',i5,2f10.3)')
     &      isp, ts, psfc
        write(h_output,'(''Obs Tb (31), Surface level calc Tb (31) '',2f10.3)')
     &      tw, ttpp(isp)
        write(h_output,'(1x)')
        do k = ltrp, isp
          write(h_output,'(i5, 4f12.4)') k, pp(k), ttpp(k), rclr_s(k), rclr_s2(k)
        enddo
        write(h_output,'(1x)')
        write(h_output,'(''Window cloud height level '',i4)') lwin
        write(h_output,'(''Window cloud fraction '',f10.3)') ecawin 
        write(h_output,'(''Window cloud height pressure '',
     &      ''and temperature '',2f10.3)') pwin,tct
c      end if
      end if

c----------------------------------------------------------------------
  
c     If using observed clear-sky radiances, make sure a valid clear-
c     sky radiance was found.

      if( (mod06ct_proc_opt .ne. 2 ) .or.
     &    (mod06ct_proc_opt .eq. 2 .and. obs_clr .eq. 1) ) then

c----------------------------------------------------------------------
  
c       Perform radiative transfer calculations for co2-slicing method.

c----------------------------------------------------------------------
  
c       Define index of first level of integration in computing RHS of
c       co2-slicing equation.
 
        is1 = isp - 1

        do k = 1,nbct

c         Compute TOA clear radiance from input profiles.

          rclr(k) = fm_modrad_emis(ts,sfc_emis(k),psfc,pp,tp,taup(1,k),
     &                             mbnds(k),isp)

c         Compute components of RHS 
          sum = 0.0
          do ll = is1,1,-1
            db = modis_planck_shift(tp(ll + 1),mbnds(k),0)
     &               - modis_planck_shift(tp(ll),mbnds(k),0)
            sum = sum - 0.5 * (taup(ll + 1,k) + taup(ll,k)) * db
            ra(ll,k) = sum

          enddo

        enddo

c----------------------------------------------------------------------
  
c       Apply clear-sky radiance bias corrections and calculate spectral
c       cloud forcing.

        num_biasflag = 0 
        do k = 1,nbct
      
c         Get bias correction for this channel
          bias(k) = cs_bias_corr(mid_px,mid_ln,mcsbnds(k))

c         Get "warm" radiance value from either calculation or observations.
          if(mod06ct_proc_opt .eq. 1) then
c           Apply clear-sky radiance bias correction.
            if((k.ne.nbct) .and. (bias(k).ne.bad_value)) then
              rtmp1 = modis_bright_shift(rclr(k),mbnds(k),0)
              rc1(k) = modis_planck_shift(rtmp1,mbnds(k),1)
              rc2(k) = rc1(k) + bias(k)
              rtmp2 = modis_bright_shift(rc2(k),mbnds(k),1)
              rwarm(k) = modis_planck_shift(rtmp2,mbnds(k),0) - rad_bias(mcsbnds(k))
              num_biasflag = num_biasflag + 1 

            else

c             No correction is available.
              rwarm(k) = rclr(k)

            end if
          else
            rwarm(k) = modis_planck(twarm(k),mbnds(k),0)
          end if

c         Get "cold" or "cloudy" radiance from observations.
          robs(k) = modis_planck(tcold(k),mbnds(k),0)

c         Calculate "cold - warm" or "clear - cloudy" value. These are 
c         components of the LHS of the co2-slicing equation.
          delr(k) = robs(k) - rwarm(k)
          spec_cloud_forcing(iout,jout,k) = delr(k)

        enddo

c----------------------------------------------------------------------
    
c       If all biases are valid, set bias_flag = 1.
        if (num_biasflag .eq. (nbct - 1)) then
          bias_flag = 1
        else
          bias_flag = 0
        end if

c----------------------------------------------------------------------
    
c       Write debug information for radiative transfer calculations.
        if(debug .gt. 0) then
c         if(line .eq. 6 .and. pixel .eq. 676) then
c         if(line .eq. 6 .and. pixel .eq. 26) then
c         if(line .eq. 6 .and. pixel .eq. 1256) then
            write(h_output,'(1x)')
            write(h_output,'(''Radiative transfer section'')')
            write(h_output,'(''Surface temperature '',f10.3)') ts
            write(h_output,'(''Values for MODIS bands '',5i10)')
     &           (mbnds(ii),ii=1,nbct)
            write(h_output,'(''Calculated TOA radiances '',5f10.3)')
     &           (rclr(ii),ii=1,nbct)
            write(h_output,'(''Calculated rc1           '',5f10.3)')
     &           (rc1(ii),ii=1,nbct)
            write(h_output,'(''Clear-sky radiance bias  '',5f10.3)')
     &           (bias(ii),ii=1,nbct)
            write(h_output,'(''Calculated rc2           '',5f10.3)')
     &           (rc2(ii),ii=1,nbct)
            write(h_output,'(''Warm or clear radiances  '',5f10.3)')
     &           (rwarm(ii),ii=1,nbct)
            write(h_output,'(''Cold or cloudy radiances '',5f10.3)')
     &           (robs(ii),ii=1,nbct)
            write(h_output,'(''Spectral cloud forcing   '',5f10.3)')
     &           (delr(ii),ii=1,nbct)
c         end if
        end if

c----------------------------------------------------------------------
  
c       Do not perform CO2-slicing for water clouds over water surfaces 
c       (check value of beta ratio).
        if(landsea .eq. 0) then
          if(beta .lt. beta_threshold) go to 3000
        endif

c----------------------------------------------------------------------
  
c       Perform CO2-slicing algorithm.

c----------------------------------------------------------------------
  
c       Compute cloud height for each channel combination (when possible).
 
        do id = 1,nsct

          krto(id) = 0
          fmsav = 1000.0
          k1 = kch(id,1)
          k2 = kch(id,2)

c         Check if values of cold minus warm are within instrument noise.
          if(delr(k1) .le. rmin(k1)) then
            if(delr(k2) .le. rmin(k2)) then
      
c               Find minimum difference between LHS and RHS of co2-slicing
c               equation. Apply emissivity correction.

                ok = .false. 
                do ll = ltrp, iw1

                  fm1 = delr(k1) / delr(k2)
                  fm2 = (ra(ll,k1) * em_adj(id)) / ra(ll,k2)
                  fm = fm1 - fm2
                  if(fm .lt. 0.0) then
                    neg = .true.
                  else
                    neg = .false.
                  end if
                  if(ll .eq. ltrp) start = neg
                  if(neg .neqv. start) then
 
                    ok = .true.
                    lev(id) = ll - 1
                    fmsav = abs(fm)
                    go to 2000
                  end if

c                 write(*,'(4i5,3f12.5,2l5)') id,k1,k2,ll,fm1,fm2,fm,neg,start

                enddo
 2000           continue

                if(fmsav .lt. 1000.0) then
                  krto(id) = 1
                end if
     
                if( (.not. ok) ) krto(id) = 0
                if( (.not. ok) .and. neg) then
                  krto(id) = 1
                  lev(id) = iw1 
                end if

            end if
          end if

        enddo

c----------------------------------------------------------------------
  
c       Write debug information for initial cloud heights.
        if(debug .gt. 0) then
          write(h_output,'(1x)')
          write(h_output,'(''Initial CO2-slicing cloud heights'')')
          write(h_output,'(''MODIS channel combinations:'')')
          write(h_output,'(8x,8(i8,''/'',i2))') (mbnds(kch(ii,1)),
     &       mbnds(kch(ii,2)),ii=1,nsct)
          write(h_output,'(''atm lev '',8i10)') (lev(ii),ii=1,nsct)
        end if

c----------------------------------------------------------------------
  
c       Compute effective cloud emissivities for every channel combination
c       for which there was a successful cloud height retrieval.
        do id = 1,nsct
          rwcld(id) = bad_value
          amo(id) = bad_value
          if(krto(id) .ne. 0) then
            ll = lev(id)
            rwcld(id) = fm_modrad_emis(tp(ll),emisrw,pp(ll),pp,tp,taup(1,kwc),
     &                             mbnds(kwc),ll)
            bot = rwcld(id) - rwarm(kwc)
            top = delr(kwc)
            if(abs(bot) .gt. 0.1) then
              ratio(id) = top / bot
              if(ratio(id) .le. 1.0 .and. ratio(id) .ge. 0.01) then
                amo(id) = ratio(id) * ecawin
              else
                krto(id) = 0
              end if
            else
              krto(id) = 0
            end if
          end if
        enddo
     
c----------------------------------------------------------------------
	  
c       Write debug information about effective emissivity calculations.
        if(debug .gt. 0) then
          write(h_output,'(1x)')
          write(h_output,'(''Initial CO2-slicing eff. emiss.'')')
          write(h_output,'(''MODIS channel combinations:'')')
          write(h_output,'(7x,8(i8,''/'',i2))') (mbnds(kch(ii,1)),
     &       mbnds(kch(ii,2)),ii=1,nsct)
          write(h_output,'(''cld rad   '',8f10.3)')
     &       (rwcld(ii),ii=1,nsct)
          write(h_output,'(''ratios    '',8f10.3)')
     &       (ratio(ii),ii=1,nsct)
          write(h_output,'(''clr rad, delr '',2f10.3)') rwarm(kwc),
     &       delr(kwc)
          write(h_output,'(''eff emiss '',8f10.3)')
     &       (amo(ii),ii=1,nsct)
        end if
  
c----------------------------------------------------------------------
  
c       Use "top-down" method for choosing a solution.

        if (platform_name(1:5) .eq. 'Aqua') then

          if(krto(1) .ne. 0) then
            cloudtop_pres_from_ratios(iout,jout,1) =
     &                          (nint(pp(lev(1)) / 5.0)) * 5.0
          end if
  
          if(krto(2) .ne. 0) then
            cloudtop_pres_from_ratios(iout,jout,2) =
     &                          (nint(pp(lev(2)) / 5.0)) * 5.0
          end if

          if(krto(4) .ne. 0) then
            cloudtop_pres_from_ratios(iout,jout,4) =
     &                          (nint(pp(lev(4)) / 5.0)) * 5.0
          end if
  
          if(krto(4) .ne. 0) then
            if(pp(lev(4)) .lt. 650.0) then
c             if(beta2 .lt. 0.90) then
                ipco2 = 4
                lco2 = lev(4)
                ngch = 1
c             end if
            end if
          end if
          if(krto(2) .ne. 0) then
            if(pp(lev(2)) .lt. 550.0) then
              ipco2 = 2
              lco2 = lev(2)
              ngch = 1
            end if
          end if
          if(krto(1) .ne. 0) then
            if(pp(lev(1)) .lt. 450.0) then
              ipco2 = 1
              lco2 = lev(1)
              ngch = 1
            end if
          end if

        else if(platform_name(1:5) .eq. 'Terra') then

          if(krto(1) .ne. 0) then
            cloudtop_pres_from_ratios(iout,jout,1) =
     &                          (nint(pp(lev(1)) / 5.0)) * 5.0
          end if
  
          if(krto(3) .ne. 0) then
            cloudtop_pres_from_ratios(iout,jout,3) =
     &                          (nint(pp(lev(3)) / 5.0)) * 5.0
          end if
  
          if(krto(3) .ne. 0) then
c           if(pp(lev(3)) .lt. 650.0) then
            if(pp(lev(3)) .lt. 650.0 .and. beta2 .lt. 0.80) then
              ipco2 = 3
              lco2 = lev(3)
              ngch = 1
            end if
          end if
          if(krto(1) .ne. 0) then
            if(pp(lev(1)) .lt. 450.0) then
              ipco2 = 1
              lco2 = lev(1)
              ngch = 1
            end if
          end if 

        end if 

c----------------------------------------------------------------------

 3000   continue

c-----------------------------------------------------------------------

c       Check if there were any valid CO2-slicing cloud top pressures.

c-----------------------------------------------------------------------

        co2_flag = 0
        if(ngch .ne. 0) then

c         Save best cloud height for output.
          pfco2 = pp(lco2)
          ecaco2 = amo(ipco2)

c         Prepare final cloud height for output and set appropritate
c         flags. Round cloud pressure to nearest 5 mb.
          pct = (nint(pfco2 / 5.0)) * 5.0
          eca = ecaco2
          ip = ipco2
          cldhgtmet_flag = ip
          tct = tp(lco2)

c         Fill output arrays.
          cloud_h_method(iout,jout) = ip * 1.0
          cloudtop_pres(iout,jout) = pct
          cloudtop_height(iout,jout) = (nint((z(lco2)*1000.0) / 50.0)) * 50.0
          cloudtop_temp(iout,jout) = tct
          cloudtop_eff_emi(iout,jout) = eca
          co2_flag = 1

        end if
  
c-----------------------------------------------------------------------

c       Perform marine low cloud height algorithm if indicated.

c-----------------------------------------------------------------------

        if ( cloudtop_pres(iout,jout) .gt. 600.0 .and. co2_flag .eq. 0) then

c         Do not perform over inland water.
          if(landsea .eq. 0 .and. (lsf .ne. 3 .and. lsf .ne. 5) ) then

            call marine_locld_retrieval(rlat, met_month, isp,
     *           ltrp, ttpp, z, tp, tw, lwin, llwin, lapse_rate)

            lco2 = llwin

c           Fill output arrays.
            pmlc = pp(llwin)
            if(pmlc .ge. psfc) pmlc = psfc
            cloudtop_pres(iout,jout) = (nint(pmlc / 5.0)) * 5.0
            cloudtop_height(iout,jout) = (nint((z(llwin)*1000.0) / 50.0)) * 50.0
            if(cloudtop_height(iout,jout) .lt. 0.0) cloudtop_height(iout,jout) = 0.0
            cloudtop_eff_emi(iout,jout) = ecawin
            cloudtop_pre_ir(iout,jout) = cloudtop_pres(iout,jout)

          end if
        end if

c----------------------------------------------------------------------
  
      end if

c----------------------------------------------------------------------
  
c     Fill in near nadir values for later output.

c----------------------------------------------------------------------
  
      if(theta .le. near_nadir_vza_limit) then
        cloudtop_height_nearnad(iout,jout) = cloudtop_height(iout,jout)
        cloudtop_pres_nearnad(iout,jout) = cloudtop_pres(iout,jout)
        cloudtop_temp_nearnad(iout,jout) = cloudtop_temp(iout,jout)
        cloudtop_eff_emi_nearnad(iout,jout) = cloudtop_eff_emi(iout,jout)
      end if

c----------------------------------------------------------------------

c     Set various flags.
  
c----------------------------------------------------------------------

c     Set quality and/or confidence flags.
      if(lco2 .eq. isp) then
        conf_flag = 1
        qual_flag = 1
      end if
      if(cloudtop_pres(iout,jout) .gt. 680.0 .and. cloudtop_eff_emi(iout,jout) .gt. 0.95
     *      .and. beta .ge. 0.95 .and. beta .lt. 900.0 .and. landsea .eq. 0) then
        conf_flag = 1
        qual_flag = 1
      end if
      if(landsea .eq. 1 .and. (psfc - 100.0) .lt. cloudtop_pres(iout,jout) ) then
        conf_flag = 1
        qual_flag = 1
      end if

c     Stratospheric cloud flag.
      if(abs(rlat) .lt. 50.0) then
        if( (tcold(2) - tcold(4)) .gt. 0.5) then
          os_top_flag = 2
        else
          os_top_flag = 1
        end if
      end if

c     High Cloud Flag.
      if(cloudtop_pres(iout,jout) .lt. 440.0) then
        hc_flag = 2
      else
        hc_flag = 1
      end if

c     Cirrus Flag.
      if(cloudtop_pres(iout,jout) .le. 680.0 .and. cloudtop_eff_emi(iout,jout) .le. 0.95) then
        ci_flag = 2
      else
        ci_flag = 1
      end if

c     QA cloud height category flag.
      if(cloudtop_pres(iout,jout) .lt. 440.0) then
        cldhgt_cat = 5
      else if(cloudtop_pres(iout,jout) .ge. 440.0 .and. cloudtop_pres(iout,jout) .lt. 680.0) then
        cldhgt_cat = 4
      else
        cldhgt_cat = 3
      end if

c----------------------------------------------------------------------
  
c     Write output debug information.

c----------------------------------------------------------------------
  
      if(debug .gt. 0) then
        write(h_output,'(1x)')
        write(h_output,'(''Output from get_co2cld'')')
        write(h_output,'(''Tropopause height '',f10.3)') 
     &    height_tropopause(iout,jout)
        write(h_output,'(''Window cloud height '',f10.3)') 
     &    cloudtop_pre_ir(iout,jout)
        write(h_output,'(''Final CO2 cloud height '',f10.3)') 
     &    cloudtop_pres(iout,jout)
        write(h_output,'(''Cloud effective emissivity '',f10.3)') 
     &    cloudtop_eff_emi(iout,jout)
        write(h_output,'(''Cloud temperature '',f10.3)') 
     &    cloudtop_temp(iout,jout)
        write(h_output,'(''CO2 cloud height solutions '',8f10.3)') 
     &    (cloudtop_pres_from_ratios(iout,jout,jj),jj=1,nsct)
        write(h_output,'(''CO2 spectral cloud forcings'',8f10.3)') 
     &    (spec_cloud_forcing(iout,jout,jj),jj=1,nsct)
        write(h_output,'(''Index of solution chosen'',f10.3)') 
     &    cloud_h_method(iout,jout)
        write(h_output,'(1x)')
        write(h_output,'(''Flag values'')')
        write(h_output,'(''qual_flag,conf_flag,co2_flag,bias_flag,ci_flag,hc_flag '',
     &    6i5)') qual_flag, conf_flag, co2_flag, bias_flag, ci_flag, hc_flag
        write(h_output,'(''os_top_flag,cldhgt_cat,nearnad_flag,cldhgtmet_flag '',
     &    4i5)') os_top_flag, cldhgt_cat, nearnad_flag, cldhgtmet_flag
      end if


c----------------------------------------------------------------------
  
      return
      end
