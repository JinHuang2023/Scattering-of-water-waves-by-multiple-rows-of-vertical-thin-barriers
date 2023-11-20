

    
    
    SUBROUTINE ASSEMB
    
    USE VAR
    IMPLICIT NONE
    
    INTEGER N,P,Q,NZ,IERR
    REAL*8 CYR(0:PTERM),CYI(0:PTERM)
    COMPLEX*16 VALUEN1,VALUEN2

    DO 10 N=-NTRUN,NTRUN
        ALPHAN(N)=WKY+2.0D0*N*PI/D
        IF(WK>=DABS(ALPHAN(N)))THEN
            GAMMAN(N)=DSQRT(WK*WK-ALPHAN(N)*ALPHAN(N))
        ELSE
            GAMMAN(N)=CI*DSQRT(ALPHAN(N)*ALPHAN(N)-WK*WK)
        END IF
10  CONTINUE

!calculate Bessel functions J_p(\alpha_n a/2)
    DO 20 N=-NTRUN,NTRUN
            
        CALL ZBESJ(ALPHAN(N)*A/2.0D0, 0.0D0, 0.0D0, 1, PTERM+1, CYR, CYI, NZ, IERR)
        BESPQN(:,N)=CYR(:)!DCMPLX(,CYI(:))
            
20  CONTINUE
    
    DO 30 Q=0,PTERM
    DO 30 P=0,PTERM
        
        N=0
        VALUEN1=1.0D0/CDSIN(GAMMAN(N)*B)*BESPQN(P,N)*BESPQN(Q,N)/GAMMAN(N)/D
        VALUEN2=1.0D0/CDTAN(GAMMAN(N)*B)*BESPQN(P,N)*BESPQN(Q,N)/GAMMAN(N)/D
        
        DO 40 N=1,NTRUN
            
            VALUEN1=VALUEN1+2.0D0*CI*CDEXP(CI*GAMMAN(N)*B)/(CDEXP(CI*2.0D0*GAMMAN(N)*B)-1.0D0)&
                *BESPQN(P,N)*BESPQN(Q,N)/GAMMAN(N)/D
            VALUEN2=VALUEN2+1.0D0/CDTAN(GAMMAN(N)*B)*BESPQN(P,N)*BESPQN(Q,N)/GAMMAN(N)/D
            
            VALUEN1=VALUEN1+2.0D0*CI*CDEXP(CI*GAMMAN(-N)*B)/(CDEXP(CI*2.0D0*GAMMAN(-N)*B)-1.0D0)&
                *BESPQN(P,-N)*BESPQN(Q,-N)/GAMMAN(-N)/D
            VALUEN2=VALUEN2+1.0D0/CDTAN(GAMMAN(-N)*B)*BESPQN(P,-N)*BESPQN(Q,-N)/GAMMAN(-N)/D
            
40      CONTINUE
        
        VALUEN2=VALUEN2+(VALT1*DCOS(WKY*A-(P+Q)*PI/2.0D0)+VALT2*DCOS((P-Q)*PI/2.0D0))

        LMAT1(Q,P)=(-CI)**DBLE(P)*CI**DBLE(Q)*VALUEN1
        LMAT2(Q,P)=(-CI)**DBLE(P)*CI**DBLE(Q)*VALUEN2

30  CONTINUE
    
    END
    
    
    