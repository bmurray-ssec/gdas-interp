      subroutine interp( nold, xold, yold, nnew, xnew, ynew )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Linearly interpolate an input array. Values outside the range of
c    the input abscissae are obtained via extrapolation.
c
c!Input Parameters:
c    NOLD    Number of elements in the input arrays
c    XOLD    Abscissa values for the input array (must be monotonic)
c    YOLD    Values of the input array
c    NNEW    Number of elements in the output arrays
c    XNEW    Abscissa values for the output array
c
c!Output Parameters:
c    YNEW    Linearly interpolated values of the output array
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
      
      integer nold, nnew
      real xold( nold ), yold( nold ), xnew( nnew ), ynew (nnew )

c ... local variables

      real slope, intercept
      integer lo, hi, j, init
      
      lo = 1
      hi = 2
      init = 1
      
      do j = 1, nnew
        
20      continue

c ...   check if output point falls between current input points
                
        if( xnew( j ) .gt. xold( hi ) ) then
          if( hi .lt. nold ) then
            lo = lo + 1
            hi = hi + 1
            init = 1
            goto 20
          endif
        endif

c ...   compute slope and intercept only when necessary

        if( init .eq. 1 ) then
          slope = ( yold( hi ) - yold( lo ) ) /
     &      ( xold( hi ) - xold( lo ) )
          intercept = yold( lo ) - slope * xold( lo )
          init = 0
        endif

c ...   compute output value

        ynew( j ) = slope * xnew( j ) + intercept

      end do
        
      end  
