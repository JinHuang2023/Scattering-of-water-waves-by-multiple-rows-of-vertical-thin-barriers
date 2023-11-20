
    MODULE VAR
    
    IMPLICIT NONE

    REAL*8,PARAMETER:: PI=4.0D0*DATAN(1.0D0),PI2=2.0D0*PI
    COMPLEX*16,PARAMETER:: CI=(0.0D0,1.0D0)

    REAL*8 THETA,WK
    REAL*8 WKX,WKY
    
    INTEGER NBAR
    REAL*8 A,B,C,D
    
    INTEGER PTERM,MTERM,NTRUN
    
    REAL*8,ALLOCATABLE:: ALPHAN(:)
    COMPLEX*16,ALLOCATABLE:: GAMMAN(:)
    COMPLEX*16,ALLOCATABLE:: BESPQN(:,:)
    
    REAL*8 VALT1,VALT2
    COMPLEX*16,ALLOCATABLE:: LMAT1(:,:),LMAT2(:,:),LMAT(:,:)
    
    COMPLEX*16,ALLOCATABLE:: BETAK(:),WTP(:,:)

    END
    
    
    
    
    
    
    PROGRAM MAIN
    
    USE VAR
    IMPLICIT NONE
    
    INTEGER M,N,NVAR
    REAL*8 WVAR1,WVAR2,DWVAR
    REAL*8 AOD
    
    OPEN(11,FILE='DATIN.TXT',      STATUS='OLD')
    
    OPEN(21,FILE='Output\OBETA.txt', STATUS='UNKNOWN')     
    OPEN(22,FILE='Output\ORT.txt',   STATUS='UNKNOWN')
        
    WRITE(21,*)'   WK        BETA'
    WRITE(22,*)'   WK        R             T'

    READ(11,*)    THETA
    READ(11,*)    WVAR1,WVAR2,DWVAR
    
    THETA=THETA/180.0D0*PI
    
    READ(11,*) A,B,C,NBAR
    READ(11,*) PTERM
    READ(11,*) MTERM
    READ(11,*) NTRUN
    
    D=A+C

    ALLOCATE (LMAT1(0:PTERM,0:PTERM),LMAT2(0:PTERM,0:PTERM))
    ALLOCATE (LMAT(0:PTERM,0:PTERM),BESPQN(0:PTERM,-NTRUN:NTRUN))
    ALLOCATE (GAMMAN(-NTRUN:NTRUN),ALPHAN(-NTRUN:NTRUN))
    ALLOCATE (BETAK(0:MTERM),WTP(0:MTERM,0:PTERM))

    !Calculate the residue of the truncation series
    AOD=2.0D0*PI*A/D
    VALT1=0.0D0
    DO N=NTRUN+1,10000000
        VALT1=VALT1+DSIN(AOD*N)/N/N
    ENDDO
    VALT2=PI*PI/6.0D0
    DO N=1,NTRUN
        VALT2=VALT2-1.0D0/N/N
    ENDDO
    
    VALT1=-VALT1*D/A/PI**3.0D0
    VALT2=-VALT2*D/A/PI**3.0D0
    
    NVAR=(WVAR2-WVAR1+1.0D-5)/DWVAR
    
    DO 100 M=0,NVAR
        WK=WVAR1+M*DWVAR
        
        WRITE(6,*) 
	    WRITE(6,*) '                   ================='
        
        WRITE(6,1111)  WK
        
        WKX=WK*DCOS(THETA)
        WKY=WK*DSIN(THETA)
        
        CALL ASSEMB
        WRITE(*,*)'   AFTER ASSEMB'

        CALL EIGENVALUE
        WRITE(*,*)'   AFTER EIGENVALUE'
        
        CALL EIGENVECTOR
        WRITE(*,*)'   AFTER EIGENVECTOR'
        
        CALL REFTRA

100 CONTINUE
  
1111 FORMAT(//,'  WAVE NUMBER=',F9.5)
     
    END 
    
    
