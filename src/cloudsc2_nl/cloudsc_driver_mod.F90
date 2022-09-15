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
  USE CLOUDSC2_ARRAY_STATE_MOD, ONLY: CLOUDSC2_ARRAY_STATE

  IMPLICIT NONE
!INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300)
INTEGER, PARAMETER :: JPRB = 8
TYPE(CLOUDSC2_ARRAY_STATE) :: GLOBAL_STATE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_PRINT( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS)
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS
  PRINT *,'NPROMA is:',NPROMA
  END SUBROUTINE CLOUDSC_DRIVER_PRINT

  SUBROUTINE CLOUDSC_DRIVER_TEST( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
     & PT, PQ, &
     & BUFFER_CML, BUFFER_LOC, &
!     & TENDENCY_CML, TENDENCY_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
!    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(10) ! cumulative tendency used for final output
!    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(10) ! local tendency from cloud scheme
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
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

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
!   if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
!   end if

    ! Global timer for the parallel region
!   CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
    !$omp& private(ZQSAT) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
!   TID = GET_THREAD_NUM()
!   CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
    print *, 'Test 1'; call flush()
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
    print *, 'Test 2'; call flush()
!         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 

    print *, 'Test 3'; call flush()
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
              &  BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
              &  BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
              &  BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL),  &
              &  BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL),  &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL))

    print *, 'Test 4'; call flush()
         ! Log number of columns processed by this thread
!         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

!      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

!     CALL TIMER%END()

!      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)
    
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

  SUBROUTINE CLOUDSC_DRIVER_INIT(NPROMA,NGPTOT,NGPTOTG)
USE PARKIND1, ONLY: JPIM
USE YOECLD   , ONLY : YRECLD
USE YOPHNC   , ONLY : YRPHNC
USE YOEPHLI  , ONLY : YREPHLI
INTEGER(KIND=JPIM) :: NGPTOT,NGPTOTG,NPROMA 
INTEGER(KIND=JPIM) ::  JK
! TODO: Create a global global memory state from serialized input data
CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)

!cloudsc2
!set up other modules not initialized in cloudsc
!YRECLD
allocate(YRECLD)
allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
! security
if (GLOBAL_STATE%KLEV>200) then
  print *, ' Dimension of ZPRES/ZPRESF is too short. '
  stop
endif
DO JK=1,GLOBAL_STATE%KLEV
  YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
ENDDO
! YRPHNC
allocate(YRPHNC)
YRPHNC%LEVAPLS2=.false.
! overload LPHYLIN
YREPHLI%LPHYLIN=.true.
  END SUBROUTINE CLOUDSC_DRIVER_INIT

END MODULE CLOUDSC_DRIVER_MOD
