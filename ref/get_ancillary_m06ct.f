      SUBROUTINE GET_ANCILLARY( LAT, LON, PRES, TEMP, MIXR, LAND,
     &  SFCTMP, PRMSL, PRSFC, PWAT, TOZONE, O3MR, UGRD, VGRD, OZONE, ICEC, SST,
     &  ICEC_SST, NISE, MET_DATE, OZN_DATE, ICE_DATE, SST_DATE,
     &  NISE_DATE, ZNBIAS_OC, ZNBIAS_DL, ZNBIAS_NL )

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c      Retrieve ancillary data items for a given latitude and longitude.
c
c !INPUT PARAMETERS:
c      LAT       Latitude (degrees, -90S to +90.0N)
c      LON       Longitude (degrees, -180W to +180E, Greenwich=0)
c
c !OUTPUT PARAMETERS:
c      PRES      Array of pressure levels (hPa)
c      TEMP      Array of atmospheric temperatures (K) at PRES(0:15)
c      MIXR      Array of water vapor mixing ratios (g/kg) at PRES(0:15)
c      LAND      Land mask (0=water, 1=land)
c      SFCTMP    Surface temperature (K)
c      PRSFC     Pressure (hPa) at surface       
c      PRMSL     Pressure (hPa) at mean sea level
c      PWAT      Precipitable water (g/cm**2)
c      TOZONE    Total column ozone (Dobsons) from NCEP
c      UGRD      Surface wind u component (m/s)
c      VGRD      Surface wind v component (m/s)
c      OZONE     TOVS Total ozone (Dobsons) **Not used**
c      ICEC      Ice concentration (fraction)
c      SST       Sea surface temperature (K) - valid over ocean only
c      ICEC_SST  Sea surface temperature ice concentration
c      NISE      NSIDC NISE snow/ice extent (see read_nise.f)
c      MET_DATE  UTC date for parameters PRES-VGRD (year,month,day,hour)
c      OZN_DATE  UTC date for parameter  OZONE     (year,month,day,hour)
c      ICE_DATE  UTC date for parameter  ICEC      (year,month,day,hour)
c      SST_DATE  UTC date for parameter  SST       (year,month,day,hour)
c      NISE_DATE UTC date for parameter  NISE      (year,month,day,hour)
c      ZNBIAS_OC Zonal mean difference between observed and modeled 
c                clear-sky radiances for ocean surfaces over 8 days
c      ZNBIAS_DL Zonal mean difference between observed and modeled 
c                clear-sky radiances for daytime land  surfaces over 8 days
c      ZNBIAS_NL Zonal mean difference between observed and modeled 
c                clear-sky radiances for nighttime land surfaces over 8 days
c    
c      The missing value for output parameters is MISSING (see below).
c
c !REVISION HISTORY:
c
c  09/12/03 R. Frey:  Added calls to pgs_pc_getreference to check for
c    the presence of grib met and sea ice files, and also for NISE files.
c    Changed error levels to "2" if any files cannot be opened or read.
c
c  12/11/02 R. Hucek: Commented code statements that access the NCEP 
c    Reynolds sst and TOVS ozone products.  
c
c  11/21/02 R. Frey:  Added subroutine read_reynsst to enable reading of
c    either formatted or unformatted Reynolds SST files.
c
c  06/04 Collection 5  R. Frey:  Added logic to acquire clear-sky
c    radiance biases.
c
c  09/05 Changed code to read zonal mean clear-sky biases instead of 25-km
c    bin data.
c
c !TEAM-UNIQUE HEADER:
c      Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !DESIGN NOTES:
c      (1) On the first call, this subroutine unpacks and reads 4 files:
c          NCEP GDAS1 meteorological analysis data,
c          NCEP TOVS total ozone data,
c          NCEP SSMI ice concentration data,
c          NCEP sea surface temperature data,
c          NSIDC NISE snow/ice extent data.
c          On subsequent calls, data is obtained from SAVEd arrays.
c
c      (2) This subroutine will not cause an exit. If errors are
c          detected (e.g. missing or bad input file), the subroutine
c          will write a 'Recoverable error' message to the LogStatus
c          file, and will return missing value(s) for the parameter(s)
c          it failed to read.
c
c      (3) No checking of data validity is done within this routine.
c          Missing data values are used only where the input file was
c          either missing or bad. The user is responsible for checking
c          that ancillary data values (e.g. SST) are within an
c          acceptable range for user's application.
c
c !END
c-----------------------------------------------------------------------
      IMPLICIT NONE

c ... Include files for PGS toolkit
      INCLUDE 'PGS_SMF.f'

c ... Include file for Ancillary product PCF numbers
      INCLUDE 'Atmos_AncData.inc'

c ... rhucek 11/17/05: Added include file mod06uw_pcfnum.inc containing 
c     MODIS product PCF numbers
      INCLUDE 'mod06uw_pcfnum.inc'
                  
c ... Input arguments
      
      REAL lat, lon

c ... Output arguments

      REAL pres(0:25), temp(0:25), mixr(0:25), o3mr(0:5),
     &  land, sfctmp, prmsl, prsfc, pwat, ugrd, vgrd, ozone, icec, sst,
     &  znbias_oc(181,5),znbias_dl(181,5),znbias_nl(181,5),tozone
      INTEGER nise, icec_sst
      INTEGER met_date(4), ozn_date(4), ice_date(4), sst_date(4),
     &  nise_date(4)
      
c ... Parameters

      REAL missing
      PARAMETER ( missing = -999.0 )
      
      INTEGER npoints_x, npoints_y
      PARAMETER ( npoints_x = 360 )
      PARAMETER ( npoints_y = 180 )

c ... Local variables

      LOGICAL INIT
      
      byte sst_ice(0:npoints_x-1,0:npoints_y-1)

      INTEGER lun, i, j, k, kk, level, ios, ret,
     &  pcfnum, reclen, status, version
         
      REAL x, x0, dx, y, y0, dy, p( 0:25 ), satmix, xlon, 
     &  met_grid( 0:359, 0:180, 0:61 ),
     &  ozn_grid( 0:359, 0:180 ),
     &  ice_grid( 0:719, 0:359 ),
     &  sst_grid( 0:npoints_x-1, 0:npoints_y-1 ), sst_bl,
     &  error_grid( 0:npoints_x-1, 0:npoints_y-1 )

      CHARACTER*3   sst_file_fmt
      CHARACTER*8   ESDT_name
      character*10 value_text
      CHARACTER*160 errmsg

      LOGICAL met_success, ozn_success, ice_success, sst_success,
     &  nise_success, csr_success
      
      CHARACTER*255 nise_file, grib_name, bias_file
      
      INTEGER gridsize
      PARAMETER ( gridsize = 721 )      
      CHARACTER*1 nise_north( gridsize, gridsize )
      CHARACTER*1 nise_south( gridsize, gridsize )

      INTEGER met_year, met_month, met_day, met_hour
      INTEGER ozn_year, ozn_month, ozn_day, ozn_hour
      INTEGER ice_year, ice_month, ice_day, ice_hour 
      INTEGER sst_year, sst_month, sst_day, sst_hour 
      INTEGER nise_year, nise_month, nise_day, nise_hour
      INTEGER nise_minute, nise_second, nise_fsecond, rtn_code

      CHARACTER*255 ECS_DateTimeGroup, HDF_AttributeName

c ... External subroutines
      EXTERNAL bl_int
      EXTERNAL blint_met
      EXTERNAL read_zonbias

c ... External functions

      INTEGER modis_grib_driver
      EXTERNAL modis_grib_driver
      
      REAL ppv
      EXTERNAL ppv
      
      INTEGER openr_temp
      EXTERNAL openr_temp
      
      INTEGER pgs_pc_getreference
      EXTERNAL pgs_pc_getreference

      INTEGER read_nise
      EXTERNAL read_nise

      INTEGER ezlh_convert
      EXTERNAL ezlh_convert
                  
c ... Save statement.
      SAVE

c ... Temperature and moisture profile pressure levels (hPa)

      DATA p / 1000.0, 975.0, 950.0, 925.0, 900.0, 850.0, 800.0,
     &   750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0,
     &   350.0, 300.0, 250.0, 200.0, 150.0, 100.0,  70.0,  50.0,
     &    30.0,  20.0,  10.0 /

c ... Initializations.

      DATA init / .true. /

c ... rhucek 11/18/05 - initialization of ancillary parameters not 
c     currently retrieved but assigned to output parameters.    
      DATA ozn_year, ozn_month, ozn_day, ozn_hour     /4*-9999/
      DATA ice_year, ice_month, ice_day, ice_hour     /4*-9999/
      DATA nise_year, nise_month, nise_day, nise_hour /4*-9999/
      DATA nise_minute, nise_second, nise_fsecond     /3*-9999/

      DATA ozn_grid /65160*-9999/
      DATA ice_grid /259200*-9999/ 

c ... rhucek 11/18/05 - character 'f' translates to decimal 102, the
c     NISE value corresponding to 'Not used'
      DATA nise_north /519841*'f'/
      DATA nise_south /519841*'f'/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     INITIALIZATION
c-----------------------------------------------------------------------

c ... Open and read input data files if this is the first call

      if ( init ) then

        do i = 1, 181
          do j = 1, 5
            znbias_oc (i,j) = -999.99
            znbias_dl (i,j) = -999.99
            znbias_nl (i,j) = -999.99
          enddo
        enddo

c ...   Set data ingest success/fail flags

        met_success = .false.
        ozn_success = .false.
        ice_success = .false.
        sst_success = .false.
        nise_success = .false.
        csr_success = .false.


c-----------------------------------------------------------------------
c      Get NCEP Meteorological Data
c-----------------------------------------------------------------------

c       Check presence of grib met file.
        version = 1
        errmsg    = ' '
        status  = pgs_pc_getreference( LUN_GDAS_0ZF, version, grib_name )
        if ( status .eq. PGS_S_SUCCESS ) then

c ...     Unpack grib met file and write to binary file
          ESDT_name = 'GDAS_0ZF'
          status = modis_grib_driver( LUN_GDAS_0ZF, LUN_TEMP_GDAS_0ZF, 
     1                                ESDT_name, errmsg,
     2                                met_year, met_month, 
     3                                met_day, met_hour)

          if ( status .ne. 0 ) then
            level = 2
            call message( 'get_ancillary', errmsg, status, level )

c ....    Open unpacked met file
          else
            pcfnum = LUN_TEMP_GDAS_0ZF
            reclen = 360*181*62*4
            status = openr_temp( pcfnum, reclen, lun )

            if ( status .ne. 0 ) then
              level = 2
              write( errmsg,'(''Error opening GDAS on PCF#'',i12)') 
     &           pcfnum
              call message( 'get_ancillary', errmsg //
     &           ' [OPERATOR ACTION: Contact SDST]', status, level )

c ......    Read the unpacked met file
            else 
              if ( status .eq. 0 ) then
                read( lun, rec = 1, iostat = ios ) met_grid

                if ( ios .ne. 0 ) then
                  level = 2
                  write( errmsg,'(''Error reading GDAS on PCF#'',i12)') 
     &               pcfnum
                  call message( 'get_ancillary', errmsg //
     &               ' [OPERATOR ACTION: Contact SDST]', ios, level )
                else
                  met_success = .true.
                endif

                call modis_io_gen_closef( pcfnum, lun )
              endif   
            endif   
          endif   

        else

          level = 1
          write( errmsg,'(''No entry found for GDAS on PCF#'',
     &       i12)') LUN_GDAS_0ZF
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )

        end if

c ... rhucek 11/18/05:  Disabled code to read SEA_ICE product
cc-----------------------------------------------------------------------
cc     Get SSM/I sea ice concentration 
cc-----------------------------------------------------------------------
cc
cc       Check presence of grib sea ice concentration file.
c        version = 1
c        ESDT_name = 'SEA_ICE'
c        errmsg    = ' '
c        status  = pgs_pc_getreference( LUN_SEA_ICE, version, grib_name )
c        if ( status .eq. PGS_S_SUCCESS ) then
c
cc ...     Unpack grib sea ice file and write to binary file
c          status = modis_grib_driver( LUN_SEA_ICE, LUN_TEMP_SEA_ICE, 
c     1                                ESDT_name, errmsg, 
c     2                                ice_year, ice_month, 
c     3                                ice_day, ice_hour ) 
c
c          if ( status .ne. 0 ) then
c            level = 2
c            call message( 'get_ancillary', errmsg, status, level ) 
c
c          else    
cc ......    Open unpacked ice file
c            pcfnum = LUN_TEMP_SEA_ICE
c            reclen = 720*360*4 
c            status = openr_temp( pcfnum, reclen, lun )
c
c            if ( status .ne. 0 ) then
c              level = 2
c              write( errmsg,'(''Error opening sea ice on PCF#'',i12)') 
c     &           pcfnum
c              call message( 'get_ancillary', errmsg //
c     &           ' [OPERATOR ACTION: Contact SDST]', status, level ) 
c
c            else    
cc ........    Read the unpacked ice file
c              if ( status .eq. 0 ) then
c                read( lun, rec = 1, iostat = ios ) ice_grid
c
c                if ( ios .ne. 0 ) then
c                  level = 2
c                  write( errmsg,'(''Error reading sea ice on PCF#'',
c     &               i12)') pcfnum
c                  call message( 'get_ancillary', errmsg //
c     &              ' [OPERATOR ACTION: Contact SDST]', ios, level ) 
c                else    
c                  ice_success = .true.
c                endif   
c
c                call modis_io_gen_closef( pcfnum, lun )
c              endif   
c            endif  
c          endif   
c
c        else
c
c          level = 1
c          write( errmsg,'(''No entry found for sea ice on PCF#'',i12)')
c     &       LUN_SEA_ICE
c          call message( 'get_ancillary', errmsg //
c     &       ' [OPERATOR ACTION: Contact SDST]', status, level )
c
c        end if


cc----------------------------------------------------------------------
cc     Get Reynolds sea surface temperature 
cc----------------------------------------------------------------------
 
         call read_reynsst ( npoints_x, npoints_y, sst_year, sst_month, 
     &                       sst_day, sst_hour, sst_grid, error_grid, 
     &                       sst_ice, sst_file_fmt, sst_success)
 
        if( (.not. sst_success) ) then
          level  =  1
          status = -1
          write( errmsg,'(''Problem in read_reynsst '')')
          call message( 'get_ancillary', errmsg //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )
        end if


c ... rhucek 11/18/05:  Disabled code to read NISE product
cc----------------------------------------------------------------------
cc     Get NISE snow extent 
cc----------------------------------------------------------------------
cc ...   Get the name of the NISE snow/ice file
c        pcfnum  = LUN_NISE
c        version = 1
c        status  = pgs_pc_getreference( pcfnum, version, nise_file )
c        if ( status .eq. PGS_S_SUCCESS ) then
c
cc ...     Read the NISE snow/ice file 
c          status = read_nise( nise_file, gridsize, 
c     &      nise_north, nise_south )
c
c          if ( status .ne. 0 ) then
c            level = 2
c            write( errmsg,'(''Error opening NISE on PCF#'',i12)') 
c     &         pcfnum
c            if ( status .lt. -1)  
c     &      write( errmsg, '(''Error reading NISE on PCF#'',i12)') 
c     &         pcfnum
c            call message( 'get_ancillary', errmsg //
c     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
c
c          else
c            version = 1
c            HDF_AttributeName = 'coremetadata.0'
c            ECS_DateTimeGroup = 'EndingDateTime'
c
c            call Parse_ECS_DateTime(pcfnum, version,
c     &                              HDF_AttributeName,
c     &                              ECS_DateTimeGroup,
c     &                              nise_year, nise_month, nise_day,
c     &                              nise_hour, nise_minute, nise_second,
c     &                              nise_fsecond, rtn_code)
c
c            if ( rtn_code .ne. 0 ) then
c              nise_success = .false.
c              level = 2
c              write( errmsg,'(''Error parsing NISE metadata on PCF#'',
c     &           i12)') pcfnum
c              call message( 'get_ancillary', errmsg //
c     &           ' [OPERATOR ACTION: Contact SDST]', status, level )
c            else
c              nise_success = .true.
c            endif
c
c          endif
c            
c        else
c
c          level = 1
c          write( errmsg,'(''No entry found for NISE on PCF#'',i12)') 
c     &       pcfnum
c          call message( 'get_ancillary', errmsg //
c     &       ' [OPERATOR ACTION: Contact SDST]', status, level )
c
c        end if

c----------------------------------------------------------------------
c     Get clear-sky radiance data. 
c----------------------------------------------------------------------

c       These data are accessed once per granule. No need to return to
c       it later.
 
c ...   rhucek 11/17/05: replaced parameter LUN_CSR8d with csb_pcfnum
c       defined in file mod06uw_pcfnum.inc; defined version = 1.
c       pcfnum  = LUN_CSR8d
        pcfnum  = csb_pcfnum
        version = 1
        status  = pgs_pc_getreference( pcfnum, version, bias_file )

        if( status .ne. PGS_S_SUCCESS ) then
          
          level = 1
          call message( 'get_ancillary',
     &      'Error getting CSR filename from PCF ' //
     &      ' [OPERATOR ACTION: Contact SDST]', status, level )

        else

c         Open and read zonal mean clear-sky radiance bias data.

          call read_zonbias(bias_file,znbias_oc,znbias_dl,znbias_nl,
     &                      csr_success)

          if( (.not. csr_success) ) then
            level = 1
            status = -1
            write( errmsg,'(''Problem in read_zonbias '')')
            call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          end if

        endif

c----------------------------------------------------------------------

c ...   Unset initialization flag
        
        init = .false.
        
      endif

c-----------------------------------------------------------------------
c     SET MISSING VALUES
c-----------------------------------------------------------------------

c ... Set return values to missing

      ozone  = missing
      icec   = missing
      sst    = missing
      nise   = int( missing )
      do i = 1, 4
        ozn_date( i )  = int( missing )
        ice_date( i )  = int( missing )
        nise_date( i ) = int( missing )
      end do
      
c-----------------------------------------------------------------------
c     GET MET AND OZONE DATA
c-----------------------------------------------------------------------

      if ( met_success .or. ozn_success ) then
      
c ...   Compute cell coordinates in met and ozn grids

        x = min( max( lon,  -179.99 ), 179.99 )
        if( x .lt. 0.0 ) x = x + 360.0
        x0 = 0.0
        dx = 1.0
        i = int( ( x - x0 + 0.5*dx ) / dx )
        if( i .eq. 360 ) i = 0
 
        y = min( max( lat, -89.99 ), 89.99 )
        y0 = 90.0
        dy = -1.0
        j = int( ( y - y0 + 0.5*dy ) / dy )

      endif
      write(*,'(''Lat, Lon, x, y: '', 2f10.3, 2i10,/)') lat, lon, i, j

      if ( met_success ) then
      
c ...     Set met return values to missing
          do k = 0, 25
            pres( k ) = missing
            temp( k ) = missing
            mixr( k ) = missing
          end do
          do k = 0, 5
            o3mr( k ) = missing
          enddo
          land   = missing
          sfctmp = missing
          prmsl  = missing
          prsfc  = missing
          pwat   = missing
          ugrd   = missing
          vgrd   = missing
          tozone   = missing
          do k = 1, 4
            met_date( k )  = int( missing )
          end do

c ...     Save output pressure levels

          do k = 0, 25
            pres( k ) = p( k )
          end do

c ...     Save output met data
c ...     (note that water vapor profile is relative humidity (%))

          do k = 0, 25
            temp( k ) = met_grid( i, j, k )
          end do
          do k = 0, 20
            mixr( k ) = met_grid( i, j, k + 26 )
          end do
          land   = met_grid( i, j, 47 )

          call blint_met( met_grid(0,0,48), i, j, y, x, sfctmp, ret )
          if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in blint_met '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          end if

          call blint_met( met_grid(0,0,49), i, j, y, x, prsfc, ret )
          if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in blint_met '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          end if
          prsfc = prsfc * 0.01

          call blint_met( met_grid(0,0,61), i, j, y, x, prmsl, ret )
          if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in blint_met '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
          end if
          prmsl = prmsl * 0.01

          pwat   = met_grid( i, j, 50 )
          ugrd   = met_grid( i, j, 51 )
          vgrd   = met_grid( i, j, 52 )

c ...     Ozone data; total column and 6 stratospheric levels. Units of 'tozone' are
c         Dobsons; for o3mr are kg/kg of dry air. Convert osmr to ppmv.
          tozone   = met_grid( i, j, 54 )
          do k = 0, 5
            kk = 55 + k
            o3mr(k) = met_grid(i, j, kk) * 1000000.0 * (28.97 / 48.0)
          enddo
          write(*,'(''GDAS pressure, temperature, RH:'')')
          do k = 0, 25
            write(*,'(3f10.3)') pres(k), temp(k), mixr(k)
          enddo
          write(*,'(1x)')
          write(*,'(''GDAS ozone, GDAS units (kg/kg) and mixing ratio (ppmv):'')')
          write(*,'(6f20.10)') (met_grid( i, j, kk ), kk=55,60),
     *       (o3mr(kk),kk=0,5)
          write(*,'(1x)')
          write(*,'(''Surface temperature and pressure:'')')
          write(*,'(2f10.3)') met_grid(0,0,48), met_grid(0,0,49)/100.0
               

c ...     Convert relative humidity profile (%) to mixing ratio (g/kg)

          do k = 0, 20
 
c ...       Compute mixing ratio at 100% relative humidity
 
            satmix = 622.0 * ppv( temp( k ) ) / pres( k )

c ...       Convert relative humidity to mixing ratio

            mixr( k ) = satmix * 0.01 * mixr( k )

          end do

c ...     Extrapolate mixing ratio profile from 100 hPa to 10 hPa

          do k = 20, 25
            mixr( k ) = max( mixr( 20 ), 0.003 ) * ( pres( k ) / 100.0 )**3
            mixr( k ) = max( mixr( k ), 0.003 )
          end do

c ...     Save date

          met_date( 1 ) = met_year
          met_date( 2 ) = met_month
          met_date( 3 ) = met_day
          met_date( 4 ) = met_hour
        
      endif
      
cc ... rhucek 12/27/05:  ozone data not used - minimize stack size by 
cc     commenting out reference to ozn_grid. 
c
c      if ( ozn_success ) then
c
cc ...   Save output ozone data
c
c        ozone  = ozn_grid( i, j )
c        
cc ...   Save date
c
c        ozn_date( 1 ) = ozn_year
c        ozn_date( 2 ) = ozn_month
c        ozn_date( 3 ) = ozn_day
c        ozn_date( 4 ) = ozn_hour
c
c      endif

c-----------------------------------------------------------------------
c     GET ICE DATA
c-----------------------------------------------------------------------

cc ... rhucek 12/27/05:  NCEP sea_ice data not used - minimize stack 
cc     size by commenting out reference to ice_grid.
c
c      if ( ice_success ) then
c      
cc ...   Compute cell coordinates in ice grid
c
c        x = min( max( lon, -179.99 ), 179.99 )
c        if( x .lt. 0.0 ) x = x + 360.0
c        x0 = 0.25
c        dx = 0.5
c        i = int( ( x - x0 + 0.5*dx ) / dx )
c        if( i .eq. 720 ) i = 0
c 
c        y = min( max( lat, -89.99 ), 89.99 )
c        y0 = 89.75
c        dy = -0.5
c        j = int( ( y - y0 + 0.5*dy ) / dy )
c
cc ...   Save output ice data
c
c        icec = ice_grid( i, j )
c
cc ...   Save date
c
c        ice_date( 1 ) = ice_year
c        ice_date( 2 ) = ice_month
c        ice_date( 3 ) = ice_day
c        ice_date( 4 ) = ice_hour
c
c      endif

c-----------------------------------------------------------------------
c     GET SST DATA
c-----------------------------------------------------------------------

      if ( sst_success ) then
                        
c ...   Compute cell coordinates in sst grid

c       Binary SST data is shifted 180 degrees relative to ASCII data.
        if(sst_file_fmt .eq. 'fmt') then
          x = min( max( lon, -179.99 ),  179.99 )
          x0 = -179.5
          xlon = lon + 180.0
        else if(sst_file_fmt .eq. 'unf') then
          if(lon .lt. 0.0) then
            xlon = lon + 360.0
          else
            xlon = lon
          end if
          x = min( max( xlon, 0.00 ),  359.99 )
          x0 = 0.5
        end if

        dx = 1.0
        i = int( ( x - x0 + 0.5*dx ) / dx )
 
        y = min( max( lat, -89.99 ), 89.99 )
        y0 = -89.5
        dy = 1.0
        j = int( ( y - y0 + 0.5*dy ) / dy )

c ...   Copy ice contration for current point into variable
        icec_sst = sst_ice(i,j)
 
c ...   Bi-linearly interpolate SST

        call bl_int(sst_grid, i, j, lat, xlon, sst_bl, ret)
        if(ret .lt. 0) then
          level = 2
          status = -1
          write( errmsg,'(''Problem in bl_int '')')
          call message( 'get_ancillary', errmsg //
     &        ' [OPERATOR ACTION: Contact SDST]', status, level )
        end if

c ...   Save output sst data

        sst = sst_bl + 273.15

c ...   Save date

        sst_date( 1 ) = sst_year
        sst_date( 2 ) = sst_month
        sst_date( 3 ) = sst_day
        sst_date( 4 ) = sst_hour

c ...   Correct Y2K problem with SST year

        if ( sst_date( 1 ) .le. 99 ) then
          if ( sst_date( 1 ) .ge. 70 ) then
            sst_date( 1 ) = sst_date( 1 ) + 1900
          else
            sst_date( 1 ) = sst_date( 1 ) + 2000
          endif
        endif
        
      endif

c-----------------------------------------------------------------------
c     GET NISE DATA
c-----------------------------------------------------------------------

cc ... rhucek 12/27/05:  NISE product not used - minimize stack size 
cc     by commenting out references to nise_south/nise_north grids.
c
c      if ( nise_success ) then
c
cc ...   Get cell coordinates for southern or northern hemisphere
cc ...   (Note the grid name strings are 'Sl' and 'Nl', L not 1)
c      
c        x = min( max( lon, -179.99 ),  179.99 )
c        y = min( max( lat, -89.99 ), 89.99 )
c
c        if ( y .lt. 0.0 ) then
c          status = ezlh_convert( 'Sl', y, x, i, j )
c        else
c          status = ezlh_convert( 'Nl', y, x, i, j )
c        endif
c
c        if ( status .ne. 0 ) then
c
c          level = 1
c          call message( 'get_ancillary', 
c     &      'Error converting lat,lon to NISE col,row' //
c     &      ' [OPERATOR ACTION: CONTACT SDST]', status, level )
c
c        else
c        
cc ...     Save output NISE data for southern or northern hemisphere
c        
c          if ( y .lt. 0.0 ) then
c            nise = ichar( nise_south( i, j ) )
c          else
c            nise = ichar( nise_north( i, j ) )
c          endif
c
c        endif
c        
cc ...   Save date
c
c        nise_date( 1 ) = nise_year
c        nise_date( 2 ) = nise_month
c        nise_date( 3 ) = nise_day
c        nise_date( 4 ) = nise_hour
c
c      endif

c-----------------------------------------------------------------------

      return
      end
