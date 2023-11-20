
    SUBROUTINE REFTRA
    
    USE VAR
    IMPLICIT NONE
    
    INTEGER P,Q,K,N,NL,NR
    REAL*8 REF,TRA
    COMPLEX*16 VALUEN(0:MTERM,0:PTERM),VALUEP1,VALUEP2,VALUEP3
    COMPLEX*16 AMS(0:MTERM,0:MTERM),AMA(0:MTERM,0:MTERM),BMS(0:MTERM),BMA(0:MTERM)
    COMPLEX*16 RN(-10:10),TN(-10:10)
    
    NL=INT(D*(WK+WKY)/2.0D0/PI+1.0D-10)
    NR=INT(D*(WK-WKY)/2.0D0/PI+1.0D-10)
    
    DO 10 P=0,PTERM
    DO 10 Q=0,MTERM
            
        N=0
        VALUEN(Q,P)=BESPQN(P,N)*BESPQN(Q,N)/GAMMAN(N)/D
        
        DO 20 N=1,NTRUN
            
            VALUEN(Q,P)=VALUEN(Q,P)+BESPQN(P, N)*BESPQN(Q, N)/GAMMAN( N)/D
            
            VALUEN(Q,P)=VALUEN(Q,P)+BESPQN(P,-N)*BESPQN(Q,-N)/GAMMAN(-N)/D
            
20      CONTINUE
        
        VALUEN(Q,P)=VALUEN(Q,P)-CI*(VALT1*DCOS(WKY*A-(P+Q)*PI/2.0D0)+VALT2*DCOS((P-Q)*PI/2.0D0))
        
        VALUEN(Q,P)=VALUEN(Q,P)*(-CI)**DBLE(P)*(-CI)**DBLE(Q)

10  CONTINUE
    
    DO 30 K=0,MTERM
        
        DO 40 Q=0,MTERM
            
            VALUEP1=(0.0D0,0.0D0)
            VALUEP2=(0.0D0,0.0D0)
            VALUEP3=(0.0D0,0.0D0)
            DO 50 P=0,PTERM
                VALUEP1=VALUEP1+WTP(K,P)*VALUEN(Q,P)
                VALUEP2=VALUEP2+WTP(K,P)*LMAT1(Q,P)
                VALUEP3=VALUEP3+WTP(K,P)*LMAT2(Q,P)
50          CONTINUE
            
            AMS(K,Q)=VALUEP1+      CDTAN(NBAR*B/2.0D0*BETAK(K))*(-1.0D0)**DBLE(Q)*(VALUEP3-CDEXP(CI*BETAK(K)*B)*VALUEP2)
            AMA(K,Q)=VALUEP1-1.0D0/CDTAN(NBAR*B/2.0D0*BETAK(K))*(-1.0D0)**DBLE(Q)*(VALUEP3-CDEXP(CI*BETAK(K)*B)*VALUEP2)
            
40      CONTINUE
        
        VALUEP1=(0.0D0,0.0D0)
        DO 60 P=0,PTERM
            VALUEP1=VALUEP1+WTP(K,P)*(-CI)**DBLE(P)*CDEXP(CI*ALPHAN(0)*A/2.0D0)*BESPQN(P,0)
60      CONTINUE
        
        BMS(K)=2.0D0*CI*VALUEP1
        BMA(K)=2.0D0*CI*VALUEP1
        
30  CONTINUE
    
!Solve the system of equations
    CALL ACGAS(AMS,BMS,MTERM+1)
    CALL ACGAS(AMA,BMA,MTERM+1)
    
    DO 70 N=-10,10
        
        RN(N)=(0.0D0,0.0D0)
        TN(N)=(0.0D0,0.0D0)
        DO 80 Q=0,MTERM
            RN(N)=RN(N)+(BMS(Q)+BMA(Q))*(-CI)**DBLE(Q)*CDEXP(-CI*ALPHAN(N)*A/2.0D0)*BESPQN(Q,N)
            TN(N)=TN(N)+(BMS(Q)-BMA(Q))*(-CI)**DBLE(Q)*CDEXP(-CI*ALPHAN(N)*A/2.0D0)*BESPQN(Q,N)
80      CONTINUE
        
        RN(N)=-RN(N)/2.0D0/CI/GAMMAN(N)/D
        TN(N)= TN(N)/2.0D0/CI/GAMMAN(N)/D
    
70  CONTINUE
    
    RN(0)=RN(0)+1.0D0
    
    REF=0.0D0
    TRA=0.0D0
    DO 100 N=-NL,NR
        REF=REF+CDABS(RN(N))**2.0D0*DREAL(GAMMAN(N))/DREAL(GAMMAN(0))
        TRA=TRA+CDABS(TN(N))**2.0D0*DREAL(GAMMAN(N))/DREAL(GAMMAN(0))
100 CONTINUE

    WRITE(22,1001)WK,REF,TRA,REF+TRA
    
1001 FORMAT(100E16.5)
     
    END