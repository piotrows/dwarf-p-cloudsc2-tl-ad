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
USE YOECLD   , ONLY : YRECLD
USE YOPHNC   , ONLY : YRPHNC
USE YOEPHLI  , ONLY : YREPHLI
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
    REAL(KIND=JPRB),    INTENT(INOUT)    :: PA(:,:,:)    ! Original Cloud fraction (t)
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
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
!         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 

print *,'PTSPHY test:',GLOBAL_STATE%PTSPHY
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
    REAL(KIND=JPRB),    INTENT(INOUT)    :: PA(:,:,:)    ! Original Cloud fraction (t)
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

  SUBROUTINE CLOUDSC_DRIVER_INIT( &
  PTSPHY, &
   NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, & 
     &  NGPTOTG, &
     & OUTPUT_PT, OUTPUT_PQ, &
     & OUTPUT_BUFFER_CML, OUTPUT_BUFFER_LOC, &
     & OUTPUT_PAP,      OUTPUT_PAPH, &
     & OUTPUT_PLU,      OUTPUT_PLUDE,    OUTPUT_PMFU,     OUTPUT_PMFD, &
     & OUTPUT_PA,       OUTPUT_PCLV,     OUTPUT_PSUPSAT,&
     & OUTPUT_PCOVPTOT )
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT,  NGPTOTG, NGPBLKS
    REAL(KIND=JPRB)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PQ(:,:,:)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCOVPTOT(:,:,:) ! Precip fraction
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSL(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSN(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSL(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSN(:,:,:)
    INTEGER(KIND=JPIM) ::  JK
    INTEGER(KIND=JPIM) ::  OUTPUT_NLEV,OUTPUT_NGPBLKS 
!Compute the number of blocks based on NGPTOT and NPROMA
OUTPUT_NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
IF(NGPBLKS.NE.OUTPUT_NGPBLKS) THEN 
   PRINT *,'Fortran NGPBLKS =',OUTPUT_NGPBLKS
   PRINT *,'Python  NGPBLKS =',NGPBLKS
   STOP 'NGPBLKS invalid in Python driver'
ENDIF
! TODO: Create a global global memory state from serialized input data
print *,'Test i1'
CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)
print *,'Test i2'
  OUTPUT_NLEV  =GLOBAL_STATE%KLEV
IF(NLEV.NE.OUTPUT_NLEV) THEN 
   PRINT *,'Fortran NLEV =',OUTPUT_NLEV
   PRINT *,'Python  NLEV =',NLEV
   STOP 'NLEV invalid in Python driver'
ENDIF

IF(PTSPHY.NE.GLOBAL_STATE%PTSPHY) THEN
print *,'PTSPHY Fortran:',GLOBAL_STATE%PTSPHY
print *,'PTSPHY PYTHON:',PTSPHY
   STOP 'PTSPHY invalid in Python driver'
ENDIF

!cloudsc2
!set up other modules not initialized in cloudsc
!YRECLD
print *,'Test i3'
allocate(YRECLD)
print *,'Test i4'
allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
! security
if (GLOBAL_STATE%KLEV>200) then
  print *, ' Dimension of ZPRES/ZPRESF is too short. '
  stop
endif
print *,'Test i5'
DO JK=1,GLOBAL_STATE%KLEV
  YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
ENDDO
print *,'Test i6'
! YRPHNC
allocate(YRPHNC)
print *,'Test i7'
YRPHNC%LEVAPLS2=.false.
print *,'Test i8'
! overload LPHYLIN
YREPHLI%LPHYLIN=.true.
print *,'Test i9'
    OUTPUT_PT(:,:,:)=          GLOBAL_STATE%PT      (:,:,:)
print *,'Test i9a'
    OUTPUT_PQ(:,:,:)=          GLOBAL_STATE%PQ      (:,:,:)
print *,'Test i9b'
    OUTPUT_BUFFER_CML(:,:,:,:)=GLOBAL_STATE%B_CML   (:,:,:,:)
print *,'Test i9c'
    OUTPUT_BUFFER_LOC(:,:,:,:)=GLOBAL_STATE%B_LOC   (:,:,:,:)
print *,'Test i9d'
      OUTPUT_PAP(:,:,:)=       GLOBAL_STATE%PAP     (:,:,:)
print *,'Test i9e'
     OUTPUT_PAPH(:,:,:)=       GLOBAL_STATE%PAPH    (:,:,:)
print *,'Test i9f'
      OUTPUT_PLU(:,:,:)=       GLOBAL_STATE%PLU     (:,:,:)
print *,'Test i9g'
    OUTPUT_PLUDE(:,:,:)=       GLOBAL_STATE%PLUDE   (:,:,:)
print *,'Test i9h'
     OUTPUT_PMFU(:,:,:)=       GLOBAL_STATE%PMFU    (:,:,:)
print *,'Test i9i'
     OUTPUT_PMFD(:,:,:)=       GLOBAL_STATE%PMFD    (:,:,:)
print *,'Test i9j'
       OUTPUT_PA(:,:,:)=       GLOBAL_STATE%PA      (:,:,:)
print *,'Test i9k'
    print *,'PCLV shape',SHAPE(GLOBAL_STATE%PCLV),SHAPE(OUTPUT_PCLV)
     OUTPUT_PCLV(:,:,:,:)=       GLOBAL_STATE%PCLV  (:,:,:,:)
print *,'Test i9l'
  OUTPUT_PSUPSAT(:,:,:)=       GLOBAL_STATE%PSUPSAT (:,:,:)
print *,'Test i9m'
 OUTPUT_PCOVPTOT(:,:,:)=       GLOBAL_STATE%PCOVPTOT(:,:,:)
print *,'Test i10'
!  OUTPUT_PFPLSL(:,:,:)=       GLOBAL_STATE%PFPLSL  (:,:,:)
!  OUTPUT_PFPLSN(:,:,:)=       GLOBAL_STATE%PFPLSN  (:,:,:)
!  OUTPUT_PFHPSL(:,:,:)=       GLOBAL_STATE%PFHPSL  (:,:,:)
!  OUTPUT_PFHPSN(:,:,:)=       GLOBAL_STATE%PFHPSN  (:,:,:)
  END SUBROUTINE CLOUDSC_DRIVER_INIT

  SUBROUTINE CLOUDSC_DRIVER_VALIDATE( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
     & INPUT_PT, INPUT_PQ, &
     & INPUT_BUFFER_CML, INPUT_BUFFER_LOC, &
     & INPUT_PAP,      INPUT_PAPH, &
     & INPUT_PLU,      INPUT_PLUDE,    INPUT_PMFU,     INPUT_PMFD, &
     & INPUT_PA,       INPUT_PCLV,     INPUT_PSUPSAT,&
     & INPUT_PCOVPTOT, &
     & INPUT_PFPLSL,   INPUT_PFPLSN,   INPUT_PFHPSL,   INPUT_PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN) :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PQ(:,:,:)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFPLSL(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFPLSN(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFHPSL(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFHPSN(:,:,:)
!   IF(INPUT_NUMOMP.NE.NUMOMP) STOP 'NUMOMP DIFFERENCE'
!   IF(INPUT_NPROMA.NE.NPROMA) STOP 'NPROMA DIFFERENCE'
!   IF(INPUT_NLEV.NE.NLEV    ) STOP 'NLEV DIFFERENCE'
!   IF(INPUT_NGPTOT.NE.NGPTOT  ) STOP 'NGPTOT DIFFERENCE'
!   IF(INPUT_NGPTOTG.NE.NGPTOTG) STOP 'NGPTOTG DIFFERENCE'
!   IF(INPUT_NGPBLKS.NE.NGPBLKS) STOP 'NGPBLKS DIFFERENCE'
!   IF(INPUT_PTSPHY.NE.PTSPHY  ) STOP 'PTSPHY  DIFFERENCE'
    print *,'PT shape',SHAPE(GLOBAL_STATE%PT),SHAPE(INPUT_PT)
    print *,'B_CML shape',SHAPE(GLOBAL_STATE%B_CML),SHAPE(INPUT_BUFFER_CML)
    print *,'PFPLSN shape',SHAPE(GLOBAL_STATE%PFPLSN),SHAPE(INPUT_PFPLSN)
    GLOBAL_STATE%PT(:,:,:)      =INPUT_PT(:,:,:)
    GLOBAL_STATE%PQ(:,:,:)      =INPUT_PQ(:,:,:)
    GLOBAL_STATE%B_CML(:,:,:,:) =INPUT_BUFFER_CML(:,:,:,:)
    GLOBAL_STATE%B_LOC(:,:,:,:) =INPUT_BUFFER_LOC(:,:,:,:)
    GLOBAL_STATE%PAP(:,:,:)     =INPUT_PAP(:,:,:)
    GLOBAL_STATE%PAPH(:,:,:)    =INPUT_PAPH(:,:,:)
    GLOBAL_STATE%PLU(:,:,:)     =INPUT_PLU(:,:,:)
    GLOBAL_STATE%PLUDE(:,:,:)   =INPUT_PLUDE(:,:,:)
    GLOBAL_STATE%PMFU(:,:,:)    =INPUT_PMFU(:,:,:)
    GLOBAL_STATE%PMFD(:,:,:)    =INPUT_PMFD(:,:,:)
    GLOBAL_STATE%PA(:,:,:)      =INPUT_PA  (:,:,:)
    GLOBAL_STATE%PCLV(:,:,:,:)    =INPUT_PCLV(:,:,:,:)
    GLOBAL_STATE%PSUPSAT(:,:,:) =INPUT_PSUPSAT(:,:,:)
    GLOBAL_STATE%PCOVPTOT(:,:,:)=INPUT_PCOVPTOT(:,:,:)
    GLOBAL_STATE%PFPLSL(:,:,:)    =INPUT_PFPLSL(:,:,:)
    GLOBAL_STATE%PFPLSN(:,:,:)    =INPUT_PFPLSN(:,:,:)
    GLOBAL_STATE%PFHPSL(:,:,:)    =INPUT_PFHPSL(:,:,:)
    GLOBAL_STATE%PFHPSN(:,:,:)    =INPUT_PFHPSN(:,:,:)
    CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)
  END SUBROUTINE CLOUDSC_DRIVER_VALIDATE

  SUBROUTINE CLOUDSC_DRIVER_FULL_FORPY( &
  PTSPHY, &
   NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, & 
     &  NGPTOTG)
USE PARKIND1, ONLY: JPIM
USE YOECLD   , ONLY : YRECLD
USE YOPHNC   , ONLY : YRPHNC
USE YOEPHLI  , ONLY : YREPHLI
    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT,  NGPTOTG, NGPBLKS
    REAL(KIND=JPRB)    :: PTSPHY       ! Physics timestep
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PT(:,:,:)    ! T at start of callpar
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PQ(:,:,:)    ! Q at start of callpar
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAP(:,:,:)   ! Pressure on full levels
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAPH(:,:,:)  ! Pressure on half levels
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLU(:,:,:)   ! Conv. condensate
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLUDE(:,:,:) ! Conv. detrained water
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFU(:,:,:)  ! Conv. mass flux up
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFD(:,:,:)  ! Conv. mass flux down
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PA(:,:,:)    ! Original Cloud fraction (t)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCLV(:,:,:,:) 
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PSUPSAT(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCOVPTOT(:,:,:) ! Precip fraction
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSL(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSN(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSL(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSN(:,:,:)
    INTEGER(KIND=JPIM) ::  JK
    INTEGER(KIND=JPIM) ::  OUTPUT_NLEV,OUTPUT_NGPBLKS 
!Compute the number of blocks based on NGPTOT and NPROMA
OUTPUT_NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
IF(NGPBLKS.NE.OUTPUT_NGPBLKS) THEN 
   PRINT *,'Fortran NGPBLKS =',OUTPUT_NGPBLKS
   PRINT *,'Python  NGPBLKS =',NGPBLKS
   STOP 'NGPBLKS invalid in Python driver'
ENDIF
! TODO: Create a global global memory state from serialized input data
CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)
  OUTPUT_NLEV  =GLOBAL_STATE%KLEV
IF(NLEV.NE.OUTPUT_NLEV) THEN 
   PRINT *,'Fortran NLEV =',OUTPUT_NLEV
   PRINT *,'Python  NLEV =',NLEV
   STOP 'NLEV invalid in Python driver'
ENDIF

IF(PTSPHY.NE.GLOBAL_STATE%PTSPHY) THEN
print *,'PTSPHY Fortran:',GLOBAL_STATE%PTSPHY
print *,'PTSPHY PYTHON:',PTSPHY
   STOP 'PTSPHY invalid in Python driver'
ENDIF

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
!!  OUTPUT_PT(:,:,:)=          GLOBAL_STATE%PT      (:,:,:)
!!  OUTPUT_PQ(:,:,:)=          GLOBAL_STATE%PQ      (:,:,:)
!!  OUTPUT_BUFFER_CML(:,:,:,:)=GLOBAL_STATE%B_CML   (:,:,:,:)
!!  OUTPUT_BUFFER_LOC(:,:,:,:)=GLOBAL_STATE%B_LOC   (:,:,:,:)
!!    OUTPUT_PAP(:,:,:)=       GLOBAL_STATE%PAP     (:,:,:)
!!   OUTPUT_PAPH(:,:,:)=       GLOBAL_STATE%PAPH    (:,:,:)
!!    OUTPUT_PLU(:,:,:)=       GLOBAL_STATE%PLU     (:,:,:)
!!  OUTPUT_PLUDE(:,:,:)=       GLOBAL_STATE%PLUDE   (:,:,:)
!!   OUTPUT_PMFU(:,:,:)=       GLOBAL_STATE%PMFU    (:,:,:)
!!   OUTPUT_PMFD(:,:,:)=       GLOBAL_STATE%PMFD    (:,:,:)
!!     OUTPUT_PA(:,:,:)=       GLOBAL_STATE%PA      (:,:,:)
!!   OUTPUT_PCLV(:,:,:,:)=       GLOBAL_STATE%PCLV  (:,:,:,:)
!!OUTPUT_PSUPSAT(:,:,:)=       GLOBAL_STATE%PSUPSAT (:,:,:)
!!UTPUT_PCOVPTOT(:,:,:)=       GLOBAL_STATE%PCOVPTOT(:,:,:)
!! OUTPUT_PFPLSL(:,:,:)=       GLOBAL_STATE%PFPLSL  (:,:,:)
!! OUTPUT_PFPLSN(:,:,:)=       GLOBAL_STATE%PFPLSN  (:,:,:)
!! OUTPUT_PFHPSL(:,:,:)=       GLOBAL_STATE%PFHPSL  (:,:,:)
!! OUTPUT_PFHPSN(:,:,:)=       GLOBAL_STATE%PFHPSN  (:,:,:)
CALL CLOUDSC_DRIVER_TEST(NUMOMP, NPROMA, GLOBAL_STATE%KLEV, NGPTOT, GLOBAL_STATE%NBLOCKS, NGPTOTG, GLOBAL_STATE%PTSPHY, &
     & GLOBAL_STATE%PT, GLOBAL_STATE%PQ, &
!    & GLOBAL_STATE%TENDENCY_CML, GLOBAL_STATE%TENDENCY_LOC, &
     & GLOBAL_STATE%B_CML,    GLOBAL_STATE%B_LOC, &
     & GLOBAL_STATE%PAP,      GLOBAL_STATE%PAPH, &
     & GLOBAL_STATE%PLU,      GLOBAL_STATE%PLUDE, &
     & GLOBAL_STATE%PMFU,     GLOBAL_STATE%PMFD, &
     & GLOBAL_STATE%PA,       GLOBAL_STATE%PCLV,  GLOBAL_STATE%PSUPSAT,&
     & GLOBAL_STATE%PCOVPTOT, &
     & GLOBAL_STATE%PFPLSL,   GLOBAL_STATE%PFPLSN,   GLOBAL_STATE%PFHPSL,   GLOBAL_STATE%PFHPSN &
     & )
  CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)
  END SUBROUTINE CLOUDSC_DRIVER_FULL_FORPY

#include "cloudsc2.F90"
#include "satur.F90"
  SUBROUTINE CLOUDSC_DRIVER_INITTEST( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
!    & PT, PQ, &
!    & BUFFER_CML, BUFFER_LOC, &
!    & PAP,      PAPH, &
!    & PLU,      PLUDE,    PMFU,     PMFD, &
!    & PA,       PCLV,     PSUPSAT,&
!    & PCOVPTOT, &
!    & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & OUTPUT_PT, OUTPUT_PQ, &
     & OUTPUT_BUFFER_CML, OUTPUT_BUFFER_LOC, &
     & OUTPUT_PAP,      OUTPUT_PAPH, &
     & OUTPUT_PLU,      OUTPUT_PLUDE,    OUTPUT_PMFU,     OUTPUT_PMFD, &
     & OUTPUT_PA,       OUTPUT_PCLV,     OUTPUT_PSUPSAT,&
     & OUTPUT_PCOVPTOT, & 
     & OUTPUT_PFPLSL,   OUTPUT_PFPLSN,   OUTPUT_PFHPSL,   OUTPUT_PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
!   REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
!   REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
!    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(10) ! cumulative tendency used for final output
!    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(10) ! local tendency from cloud scheme
!   REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
!   REAL(KIND=JPRB), INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
!   REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
!   REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
!   REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
!   REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
!   REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
!   REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
!   REAL(KIND=JPRB),    INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
!   REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
!   REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
!   REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
!   REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
!   REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
!   REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
!   REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PQ(:,:,:)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSL(:,:,:)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSN(:,:,:)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSL(:,:,:)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSN(:,:,:)

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND
    INTEGER(KIND=JPIM) ::  JK

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

!Compute the number of blocks based on NGPTOT and NPROMA
!OUTPUT_NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
!IF(NGPBLKS.NE.OUTPUT_NGPBLKS) THEN 
!   PRINT *,'Fortran NGPBLKS =',OUTPUT_NGPBLKS
!   PRINT *,'Python  NGPBLKS =',NGPBLKS
!   STOP 'NGPBLKS invalid in Python driver'
!ENDIF
!! TODO: Create a global global memory state from serialized input data
 CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)
!  OUTPUT_NLEV  =GLOBAL_STATE%KLEV
!IF(NLEV.NE.OUTPUT_NLEV) THEN 
!   PRINT *,'Fortran NLEV =',OUTPUT_NLEV
!   PRINT *,'Python  NLEV =',NLEV
!   STOP 'NLEV invalid in Python driver'
!ENDIF

IF(PTSPHY.NE.GLOBAL_STATE%PTSPHY) THEN
!print *,'PTSPHY Fortran:',GLOBAL_STATE%PTSPHY
!print *,'PTSPHY PYTHON:',PTSPHY
!   STOP 'PTSPHY invalid in Python driver'
ENDIF

print *,'Test a'
!cloudsc2
!set up other modules not initialized in cloudsc
!YRECLD
allocate(YRECLD)
allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
print *,'Test b'
! security
if (GLOBAL_STATE%KLEV>200) then
  print *, ' Dimension of ZPRES/ZPRESF is too short. '
  stop
endif
DO JK=1,GLOBAL_STATE%KLEV
  YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
ENDDO
print *,'Test c'
! YRPHNC
allocate(YRPHNC)
YRPHNC%LEVAPLS2=.false.
print *,'Test d'
! overload LPHYLIN
YREPHLI%LPHYLIN=.true.
1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
!   if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
!   end if
print *,'Test e'
CALL CLOUDSC_DRIVER_TEST(NUMOMP, NPROMA, GLOBAL_STATE%KLEV, NGPTOT, GLOBAL_STATE%NBLOCKS, NGPTOTG, GLOBAL_STATE%PTSPHY, &
     & GLOBAL_STATE%PT, GLOBAL_STATE%PQ, &
!    & GLOBAL_STATE%TENDENCY_CML, GLOBAL_STATE%TENDENCY_LOC, &
     & GLOBAL_STATE%B_CML,    GLOBAL_STATE%B_LOC, &
     & GLOBAL_STATE%PAP,      GLOBAL_STATE%PAPH, &
     & GLOBAL_STATE%PLU,      GLOBAL_STATE%PLUDE, &
     & GLOBAL_STATE%PMFU,     GLOBAL_STATE%PMFD, &
     & GLOBAL_STATE%PA,       GLOBAL_STATE%PCLV,  GLOBAL_STATE%PSUPSAT,&
     & GLOBAL_STATE%PCOVPTOT, &
     & GLOBAL_STATE%PFPLSL,   GLOBAL_STATE%PFPLSN,   GLOBAL_STATE%PFHPSL,   GLOBAL_STATE%PFHPSN &
     & )
print *,'Test f'
!    ! Global timer for the parallel region
!!   CALL TIMER%START(NUMOMP)
!
!    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
!    !$omp& private(ZQSAT) &
!    !$omp& num_threads(NUMOMP)
!
!    ! Local timer for each thread
!!   TID = GET_THREAD_NUM()
!!   CALL TIMER%THREAD_START(TID)
!
!    !$omp do schedule(runtime)
!    DO JKGLO=1,NGPTOT,NPROMA
!       IBL=(JKGLO-1)/NPROMA+1
!       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
!
!         !-- These were uninitialized : meaningful only when we compare error differences
!         PCOVPTOT(:,:,IBL) = 0.0_JPRB
!!         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB
!
!         ! Fill in ZQSAT
!         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
!              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 
!
!         CALL CLOUDSC2 ( &
!              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
!              & PTSPHY,&
!              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
!              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
!              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
!              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
!!             &  TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
!!             &  TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
!!             &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
!!             &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
!              &  BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
!              &  BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
!              &  BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL),  &
!              &  BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL),  &
!              &  PSUPSAT(:,:,IBL), &
!              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
!              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL))
!
!         ! Log number of columns processed by this thread
!!         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
!      ENDDO
!
!      !-- The "nowait" is here to get correct local timings (tloc) per thread
!      !   i.e. we should not wait for slowest thread to finish before measuring tloc
!      !$omp end do nowait
!
!!      CALL TIMER%THREAD_END(TID)
!
!      !$omp end parallel
!
!!     CALL TIMER%END()
!
!!      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)
print *,'Test g'
    OUTPUT_PT(:,:,:)=          GLOBAL_STATE%PT      (:,:,:)
    OUTPUT_PQ(:,:,:)=          GLOBAL_STATE%PQ      (:,:,:)
    OUTPUT_BUFFER_CML(:,:,:,:)=GLOBAL_STATE%B_CML   (:,:,:,:)
    OUTPUT_BUFFER_LOC(:,:,:,:)=GLOBAL_STATE%B_LOC   (:,:,:,:)
      OUTPUT_PAP(:,:,:)=       GLOBAL_STATE%PAP     (:,:,:)
     OUTPUT_PAPH(:,:,:)=       GLOBAL_STATE%PAPH    (:,:,:)
      OUTPUT_PLU(:,:,:)=       GLOBAL_STATE%PLU     (:,:,:)
    OUTPUT_PLUDE(:,:,:)=       GLOBAL_STATE%PLUDE   (:,:,:)
     OUTPUT_PMFU(:,:,:)=       GLOBAL_STATE%PMFU    (:,:,:)
     OUTPUT_PMFD(:,:,:)=       GLOBAL_STATE%PMFD    (:,:,:)
       OUTPUT_PA(:,:,:)=       GLOBAL_STATE%PA      (:,:,:)
     OUTPUT_PCLV(:,:,:,:)=       GLOBAL_STATE%PCLV  (:,:,:,:)
  OUTPUT_PSUPSAT(:,:,:)=       GLOBAL_STATE%PSUPSAT (:,:,:)
 OUTPUT_PCOVPTOT(:,:,:)=       GLOBAL_STATE%PCOVPTOT(:,:,:)
   OUTPUT_PFPLSL(:,:,:)=       GLOBAL_STATE%PFPLSL  (:,:,:)
   OUTPUT_PFPLSN(:,:,:)=       GLOBAL_STATE%PFPLSN  (:,:,:)
   OUTPUT_PFHPSL(:,:,:)=       GLOBAL_STATE%PFHPSL  (:,:,:)
   OUTPUT_PFHPSN(:,:,:)=       GLOBAL_STATE%PFHPSN  (:,:,:)
print *,'Test h'
    
  CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)
  END SUBROUTINE CLOUDSC_DRIVER_INITTEST
END MODULE CLOUDSC_DRIVER_MOD
!
