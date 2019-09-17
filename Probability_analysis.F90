PROGRAM WSMTD_rw_2D
IMPLICIT NONE
REAL*8 ::  den,alpha,Ro,kt0,kt,ktb,dum,s1
REAL*8, ALLOCATABLE :: cv1(:),prob(:),fes(:),fes1(:),grid(:)
REAL*8, ALLOCATABLE :: dummy(:,:),cv(:,:)
REAL*8, ALLOCATABLE :: gridmin(:),gridmax(:),griddif(:)
INTEGER,ALLOCATABLE :: nbin(:)
INTEGER :: md_steps,dummy1,i,j,index1,k,t_min,t_max,nbin1,i_md,i_s1,narg,w_cv
INTEGER :: ii,ncv
LOGICAL :: pmf,inpgrid
CHARACTER*120 :: arg 
REAL*8, PARAMETER :: kb = 1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: kj_to_kcal = 0.239006 

OPEN(11,FILE='cvmdck_mtd',STATUS='unknown')
OPEN(12,FILE='cv.dat',STATUS='unknown')
!      OPEN(14,FILE='cvmdck_mtd',STATUS='unknown')

CALL get_steps(11,md_steps)
!

kt0 = 300.D0  ; kt = 300.D0
t_min = 1     ; t_max = md_steps
pmf = .FALSE. ; inpgrid = .false.
t_max = md_steps ; narg = IARGC()

DO i=1,narg
  CALL GETARG(i,arg)
  IF(INDEX(arg,'-T0').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt0
  ELSEIF(INDEX(arg,'-T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt
  ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_min
  ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_max
     IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
   ELSEIF(INDEX(arg, '-ncv') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)ncv
  ELSE IF(INDEX(arg,'-grid').NE.0)THEN

ALLOCATE(gridmin(ncv),gridmax(ncv),griddif(ncv))
j = 0 
DO ii = 1,ncv
j = j + 1
      CALL GETARG(i+j,arg)
      READ(arg,*)gridmin(ii)
      CALL GETARG(i+j+1,arg)
      READ(arg,*)gridmax(ii)
      CALL GETARG(i+j+2,arg)
      READ(arg,*)griddif(ii)
j = j + 2 
WRITE(*,*)ii,gridmin(ii),gridmax(ii),griddif(ii)
ENDDO
      inpgrid=.true.
  ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_cv
  END IF
END DO

!      md_steps=t_max
   WRITE(*,'(A,I10)')'No: of MD  steps            =',md_steps
   WRITE(*,'(A,I10)')'No: of max MD  steps        =',t_max
   WRITE(*,'(A,I10)')'No: of min MD  steps        =',t_min
   WRITE(*,'(A,I10)')'No: of CV                   =',ncv

ALLOCATE(cv(ncv,md_steps))
ALLOCATE(dummy(ncv,md_steps))
ALLOCATE(nbin(ncv))

101 FORMAT (I5,4F16.6)
DO i_md=1,md_steps
   READ(11,*)dummy1,dummy1,(dummy(j,i_md), j=1,ncv),(cv(j,i_md) ,j=1,ncv)
!   WRITE(*,*)dummy1,(dummy(j,i_md), j=1,ncv),(cv(j,i_md) ,j=1,ncv)
   WRITE(12,101)dummy1,(cv(j,i_md) ,j=1,ncv)
END DO

   nbin(1) = NINT((gridmax(1)-gridmin(1))/griddif(1)) + 1
   WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
   WRITE(*,'(A10,3F8.4,I10)')'US  COORD:', gridmin(1),gridmax(1),griddif(1),nbin(1)

ALLOCATE(prob((nbin(1))))

!calculate prob        
  den=0.d0 ; prob=0.d0
 DO i_md=1,md_steps
    IF ((i_md.ge.t_min).and.(i_md.le.t_max)) THEN
     index1 = nint((cv(1,i_md)-gridmin(1))/griddif(1)) +1
     prob(index1) = prob(index1) + 1.d0
    END IF
  END DO
  
  DO index1=1,nbin(1)
     den=den+prob(index1)
  END DO   
OPEN(2,FILE='Pu.dat',STATUS='unknown')
 DO i_s1=1,nbin(1)
      s1=DFLOAT(i_s1-1)*griddif(1)+gridmin(1)
      WRITE(2,*)s1,1.0,prob(i_s1)/(den*griddif(1))
!      WRITE(*,*)s1,1.0,prob(i_s1)/(den*griddif(1))
 END DO
WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'

END PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
SUBROUTINE get_steps(iunit,nsteps)
IMPLICIT NONE
INTEGER iunit, nsteps
INTEGER ios
nsteps=0
REWIND(iunit)
Read_Loop: DO
   READ(iunit,*,IOSTAT=ios)
   IF(ios.ne.0)EXIT Read_Loop
   nsteps=nsteps+1
END DO Read_Loop 
REWIND(iunit)
END 
!---------------------------------------------------------------------!
SUBROUTINE check_files(iunit,dt)
IMPLICIT NONE
INTEGER iunit, dt
INTEGER ii, jj,i,ios
dt = 0 ;  i = 2
REWIND(iunit)
READ(iunit,*)ii
READ(iunit,*)jj
dt=jj-ii
ii=jj
RLoop: DO 
  i=i+1
  READ(iunit,*,IOSTAT=ios)jj
  IF(ios.ne.0)EXIT RLoop
  IF(jj.ne.ii+dt)THEN
     print *, '!!ERROR: Steps are not at constant stride!!'
     print *, '!!       Unit No:',iunit,'!!'
     print *, '!!       Line No:',i,'!!'
     print *, '!! Expected stride =', dt,'!!'
     print *, '!! Actual stride =', jj-ii,'!!'
     STOP
  END IF
  ii=jj
END DO RLoop
REWIND(iunit)
END 
!---------------------------------------------------------------------!
SUBROUTINE get_gridmin_max(iunit,gridmin1,gridmax1,griddif1,gridmin2,gridmax2,griddif2)
IMPLICIT NONE 
INTEGER :: iunit 
REAL*8  :: gridmin1, gridmax1, griddif1 ,gridmin2, gridmax2, griddif2
INTEGER :: ii, ios
REAL*8  :: cv1, cv2
INTEGER, PARAMETER :: Def_Grid_Size=101
REWIND(iunit)
READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
if(ios.ne.0)stop 'ERROR reading CV.dat'
gridmin1 = cv1 ; gridmax1 = cv1
gridmin2 = cv2 ; gridmax2 = cv2
RLoop: DO 
  READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
  if(ios.ne.0)EXIT RLoop
  gridmin1=MIN(gridmin1,cv1) ;  gridmin2=MIN(gridmin2,cv2)
  gridmax1=MAX(gridmax1,cv1) ;  gridmax2=MAX(gridmax2,cv2)
END DO RLoop
griddif1=(gridmax1-gridmin1)/DFLOAT(Def_Grid_Size)
griddif2=(gridmax2-gridmin2)/DFLOAT(Def_Grid_Size)
END
!---------------------------------------------------------------------!

