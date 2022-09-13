! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CLOUDSC_DRIVER_MOD
!  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, NCLDQL, NCLDQI 
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM
  USE EC_PMON_MOD, ONLY: EC_PMON

  IMPLICIT NONE
!INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300)
INTEGER, PARAMETER :: JPRB = 8
CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_PRINT( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS)
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS
  PRINT *,'NPROMA is:',NPROMA
  END SUBROUTINE CLOUDSC_DRIVER_PRINT

  SUBROUTINE CLOUDSC_DRIVER_TEST( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS, PTSPHY, &
     & PT, PQ, &
!    & TENDENCY_CML, TENDENCY_LOC, &
     & GLOBAL_STATE_CML_u  , GLOBAL_STATE_CML_v, GLOBAL_STATE_CML_T, &
     & GLOBAL_STATE_CML_o3 , GLOBAL_STATE_CML_q, GLOBAL_STATE_CML_a, &
     & GLOBAL_STATE_CML_cld, &
     & GLOBAL_STATE_LOC_u  , GLOBAL_STATE_LOC_v, GLOBAL_STATE_LOC_T, &
     & GLOBAL_STATE_LOC_o3 , GLOBAL_STATE_LOC_q, GLOBAL_STATE_LOC_a, &
     & GLOBAL_STATE_LOC_cld, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
!   TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(NBLOCKS) ! cumulative tendency used for final output
!   TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(NBLOCKS) ! local tendency from cloud scheme
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_u (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_v (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_T (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_o3(:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_q (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_a (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(IN)    :: TENDENCY_CML_cld (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_u (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_v (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_T (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_o3(:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_q (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_a (:,:,:,:)    ! 
!   REAL(KIND=JPRB),    INTENT(OUT)   :: TENDENCY_LOC_cld (:,:,:,:)    ! 
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

    REAL(KIND=JPRB)    :: LOC_u(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: LOC_v(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: LOC_T(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: LOC_o3(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: LOC_q(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: LOC_a(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: LOC_cld(NPROMA,NLEV,NCLV) ! local array
    REAL(KIND=JPRB)    :: CML_u(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: CML_v(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: CML_T(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: CML_o3(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: CML_q(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: CML_a(NPROMA,NLEV) ! local array
    REAL(KIND=JPRB)    :: CML_cld(NPROMA,NLEV,NCLV) ! local array

    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_u  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_v  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_T  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_o3 (NPROMA,NLEV,NBLOCKS)
    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_q  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_a  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_LOC_cld(NPROMA,NLEV,NCLV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_u  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_v  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_T  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_o3 (NPROMA,NLEV,NBLOCKS)
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_q  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_a  (NPROMA,NLEV,NBLOCKS) 
    REAL(KIND=JPRB)    :: GLOBAL_STATE_CML_cld(NPROMA,NLEV,NCLV,NBLOCKS) 
    
    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
    !$omp& private(ZQSAT) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
!ZP        TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB
         LOC_cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 
!ZP         LOC_T(:,:)=TENDENCY_LOC(IBL)%T
!        LOC_T(:,:)=GLOBAL_STATE_LOC_T(:,:,IBL)
!        LOC_Q(:,:)=GLOBAL_STATE_LOC_Q(:,:,IBL)
         CML_T(:,:)=GLOBAL_STATE_CML_T(:,:,IBL)
         CML_Q(:,:)=GLOBAL_STATE_CML_Q(:,:,IBL)
         CML_CLD(:,:,:)=GLOBAL_STATE_CML_CLD(:,:,:,IBL)

         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY,&
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
!             &  TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
!             &  TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
!             &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
!             &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
              &  LOC_T  (:,:)       , CML_T  (:,:)       , LOC_Q  (:,:)       , CML_Q  (:,:)       , &
              &  LOC_CLD(:,:,NCLDQL), CML_CLD(:,:,NCLDQL), LOC_CLD(:,:,NCLDQI), CML_CLD(:,:,NCLDQI), &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL))

         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)
    
  END SUBROUTINE CLOUDSC_DRIVER_TEST

  SUBROUTINE CLOUDSC_DRIVER( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS, PTSPHY, &
     & PT, PQ, &
     & TENDENCY_CML, TENDENCY_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(10) ! cumulative tendency used for final output
    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(10) ! local tendency from cloud scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
    !$omp& private(ZQSAT) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 

         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY,&
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
              &  TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
              &  TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
              &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
              &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL))

         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)
    
  END SUBROUTINE CLOUDSC_DRIVER

END MODULE CLOUDSC_DRIVER_MOD
