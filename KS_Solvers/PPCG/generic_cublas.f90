!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_threaded_memset(array, val, length)
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  REAL(DP), INTENT(OUT) :: array(length)
#if defined(__CUDA)
  attributes(device) :: array
#endif
  REAL(DP), INTENT(IN) :: val
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
!$cuf kernel do(1)  
  DO i=1, length
     array(i) = val
  ENDDO
  !
END SUBROUTINE gpu_threaded_memset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_threaded_assign(array_out, array_in, kdimx, nact, use_idx, idx, bgrp_root_only)
  !  assign (copy) a complex array in a threaded way
  !  array_out( 1:kdimx, 1:nact ) = array_in( 1:kdimx, 1:nact )       or
  !  array_out( 1:kdimx, 1:nact ) = array_in( 1:kdimx, idx(1:nact) )
  !  if the index array idx is given
  !  if  bgrp_root_only is present and .true. the assignement is made only by the 
  !  MPI root process of the bgrp and array_out is zeroed otherwise
#if defined(__CUDA)
  USE cudafor
#endif
  USE util_param,   ONLY : DP
  USE mp_bands_util,      ONLY : root_bgrp_id, nbgrp, my_bgrp_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: kdimx, nact
  COMPLEX(DP), INTENT(OUT) :: array_out( kdimx, nact )
  COMPLEX(DP), INTENT(IN)  :: array_in ( kdimx, * )
  INTEGER, INTENT(IN) :: idx( * )
#if defined(__CUDA)
  attributes(device) :: array_out, array_in, idx
#endif
  LOGICAL, INTENT(IN) :: bgrp_root_only
  LOGICAL, INTENT(IN) :: use_idx
  !
  INTEGER, PARAMETER :: blocksz = 256
  INTEGER :: nblock
  
  INTEGER :: i, j
  !
  IF (kdimx <=0 .OR. nact<= 0) RETURN
  !
  IF (bgrp_root_only .AND. ( my_bgrp_id /= root_bgrp_id ) ) THEN
     call threaded_memset( array_out, 0.d0, 2*kdimx*nact )
     RETURN
  END IF
  
  nblock = (kdimx - 1)/blocksz  + 1
  
  IF (use_idx ) THEN
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) 
       ENDDO 
     ENDDO
  ELSE
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) 
       ENDDO 
     ENDDO
  END IF
  !
END SUBROUTINE gpu_threaded_assign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_threaded_backassign(array_out, idx, array_in, kdimx, nact, use_a2, a2_in )
  !  assign (copy) a complex array in a threaded way
  !  array_out( 1:kdimx, idx(1:nact) ) = array_in( 1:kdimx, 1:nact )      or
  !  array_out( 1:kdimx, idx(1:nact) ) = array_in( 1:kdimx, 1:nact )  + a2_in( 1:kdimx, idx(1:nact) (
  !  if a2_in is present
  !  the index array idx is mandatory otherwise one could use previous routine)
#if defined(__CUDA)
  USE cudafor
#endif
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: kdimx, nact
  COMPLEX(DP), INTENT(INOUT) :: array_out( kdimx, * ) ! we don't want to mess with un referenced columns
  COMPLEX(DP), INTENT(IN)  :: array_in ( kdimx, nact )
  COMPLEX(DP), INTENT(IN)  :: a2_in ( kdimx, * )
  INTEGER, INTENT(IN)      :: idx( * )
#if defined(__CUDA)
  attributes(device) :: array_out, array_in, a2_in, idx
#endif
  LOGICAL, INTENT(IN) :: use_a2
  !
  INTEGER, PARAMETER :: blocksz = 256
  INTEGER :: nblock

  INTEGER :: i, j
  !
  IF (kdimx <=0 .OR. nact<= 0) RETURN
  !

  nblock = (kdimx - 1)/blocksz  + 1

  IF ( use_a2) THEN
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) = &
            array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i )     + & 
                a2_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) 
       ENDDO 
     ENDDO
  ELSE
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) 
       ENDDO 
     ENDDO
  END IF
  !
END SUBROUTINE gpu_threaded_backassign
