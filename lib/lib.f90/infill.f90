!-----------------------------------------------------------------------------
!
! infill.f90 -- Fill gaps in sequences by simple interpolation.
!
! 2023-sep-28	Original version.  By Dave Allured, NOAA/PSL/CIRES.
! 2023-sep-29	Add diagnostic verbosity control and ID labels.
!
! This version locates gaps (missing values) in sequences, and
! fills them by linear interpolation.  User specifies the maximum
! gap size to fill.  Larger gaps are not filled.
!
! Only complete gaps with data on both sides are infilled.
! Terminal gaps on either end are not infilled.
!
! This version operates on multiple sequences in a 2-D array (t,x).
! Infill is performed on the left hand dimension (t), which is the
! contiguous or fastest moving dimension in fortran memory order.
! This is crafted so as to auto-parallelize easily, when
! appropriate.
!
! For efficiency, this version modifies the data array in place,
! rather than making a copy.
!
!-----------------------------------------------------------------------------

module infill_mod
contains

subroutine infill (x, vmiss, max_length, diag, labels)
   implicit none

   integer, parameter :: dp = kind (1d0)	! kind = double precision

   real(dp),     intent(inout)        :: x(:,:)	    ! 2-D input array (t,k)
   real(dp),     intent(in)           :: vmiss	    ! missing value code
   integer,      intent(in)           :: max_length ! maximum gap length to fill
   						    ! (number of missing values)
   integer,      intent(in)           :: diag	    ! diagnostic verbosity
   character(*), intent(in), optional :: labels(:)  ! ID labels for K dimension

! Local variables.

   character label2*99

   integer i, istart, ilast, nflips, ngaps, n_elems
   integer k, t, t1, t2, nk, nt, llen

   real(dp) full_step, stride

   integer, allocatable :: inds(:)

   logical, allocatable :: missing(:)
   logical, allocatable :: flip(:)

! Main loop over each input sequence.

   nt = size (x, 1)					! get dimensions
   nk = size (x, 2)

k_loop: &
   do k = 1, nk

! Setup for optional diagnostic labels.

      if (diag >= 3) then
         if (present (labels)) then
            label2 = labels(k)			! use caller's labels
            llen   = len (labels)		! match caller's string length
         else
            write (label2, '(i4)') k		! if labels not provided...
            llen   = 4				! print index numbers instead
         end if
      end if

! Get position of each transition; non-missing to missing, and reverse.

      missing = (x(:,k) .eq. vmiss)			! boolean array (t)

      flip = (missing(2:nt) .neqv. missing(1:nt-1))	! transition mask
      							! boolean, size (t-1)

      inds = pack ((/ (t, t = 2, nt) /), flip(:))	! positions of all flips

! Check and skip this sequence, if no transitions.

      nflips = size (inds)		! can be size zero

      if (diag >= 4 .and. nflips == 0) print '(2a,99i5)', label2(1:llen), &
         '  nmiss, fcount, nflips =', count (missing), count (flip), nflips

      if (nflips == 0) cycle k_loop	! skip if no missing or all missing

! Find the first first transition from non-missing to missing.

      if (x(inds(1),k) == vmiss) then
         istart = 1			! sequence starts with data present
      else
         istart = 2			! sequence starts with missing values
      end if

! Compute number of complete gaps (present -- missing -- present).
! Integer division truncates, ignores any final, unterminated block.

      ngaps = (nflips - istart + 1) / 2

! Loop over all gaps in the current sequence.
! Stride by two through the alternating flip indices.

      ilast = istart + (2 * (ngaps - 1))

      if (diag >= 4) print '(2a,99i5)', label2(1:llen), &
         '  nmiss, fcount, nflips, ngaps, istart, ilast =', &
         count(missing), count(flip), nflips, ngaps, istart, ilast

gap_loop: &
      do i = istart, ilast, 2
         t1      = inds(i)		! first element position in gap
         t2      = inds(i+1) - 1	! last  element position in gap
         n_elems = (t2 - t1) + 1	! number of elements (gap size)

         if (diag >= 4 .and. n_elems > max_length) print '(2a,99i5)', &
            label2(1:llen), '  i, t1, t2, n_elems =', t1, t2, n_elems

         if (n_elems > max_length) cycle gap_loop   ! skip gap if too large

! Infill by linear interpolation through the current gap.
! Compute from the two adjacent values OUTSIDE the gap.
! Overwrite missing values.

         full_step = x(t2+1,k) - x(t1-1,k)	! value change over gap, signed
         stride    = full_step / (n_elems + 1)	! value change at each element

         if (diag >= 3) print *

         if (diag >= 4) print '(2a,4i5,2f9.2)', label2(1:llen), &
            '  i, t1, t2, n_elems, stride, full_step =', &
               i, t1, t2, n_elems, stride, full_step

         if (diag >= 3) print '(a,4i5,f9.2,a,99f9.2)', label2(1:llen), &
               i, t1, t2, n_elems, stride, '  before:', x(t1-1:t2+1,k)

         x(t1:t2,k)  = x(t1-1,k) + (stride * (/ (t, t = 1, n_elems) /))

         if (diag >= 3) print '(a,4i5,f9.2,a,99f9.2)', label2(1:llen), &
               i, t1, t2, n_elems, stride, '  before:', x(t1-1:t2+1,k)

      end do gap_loop

   end do k_loop

! All sequences completed.  Return with input array modified.

end subroutine infill
end module infill_mod
