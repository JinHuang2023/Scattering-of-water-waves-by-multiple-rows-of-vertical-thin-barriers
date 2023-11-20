
    
    
    
    
    SUBROUTINE EIGENVALUE
    
    USE VAR
    IMPLICIT NONE
    
    INTEGER M,P,IERR
    real*8 ar(0:PTERM,0:PTERM),ai(0:PTERM,0:PTERM)
    real*8 fv1(0:PTERM),fv2(0:PTERM),fv3(0:PTERM)
    real*8 wr(0:PTERM),wi(0:PTERM),zr(0:PTERM,0:PTERM),zi(0:PTERM,0:PTERM)
    COMPLEX*16 LM(0:PTERM,0:PTERM),BK(0:MIN(MTERM+2,PTERM))

    !Calculate the inverse of LMAT1
    LM=LMAT2
    CALL INVERSE(LM(0:pterm,0:pterm),LMAT,PTERM+1,PTERM+1)

    LMAT=MATMUL(LMAT1(0:pterm,0:pterm),LMAT)
    
    !Calaulate the eigenvalue and eigenvector of LMAT
    ar=dreal(LMAT)
    ai=dimag(LMAT)
    CALL cgK(PTERM+1,PTERM+1,ar,ai,wr,wi,0,zr,zi,fv1,fv2,fv3,ierr)
    
    wr(:)=1.0d0/wr(:)
    
    !Sort the eigenvalue
    CALL PIKER_R(PTERM+1,WR)
    
    M=-1
    DO 10 P=0,PTERM
        IF(WR(P)<-1.0D0)THEN
            M=M+1
            BK(M)=(PI+CI*DACOSH(-WR(P)))/B
        ELSE IF(WR(P)<=1.0D0)THEN
            M=M+1
            BK(M)=DACOS(WR(P))/B
        ELSE
            M=M+1
            BK(M)=CI*DACOSH(WR(P))/B
        END IF
        IF(M==MTERM+2)GOTO 20
10  CONTINUE

20  CONTINUE
    
    !Sort the imaginary part of BK
    CALL PIKER(MIN(MTERM+3,PTERM+1),BK)
    
    BETAK(0:MTERM)=BK(0:MTERM)
    
    WRITE(21,1001)WK,(BETAK(M)*B,M=0,MTERM)
    
1001 FORMAT(F16.5,20F16.5)  
     
    END

    
    
    
    SUBROUTINE EIGENVECTOR
    
    USE VAR
    IMPLICIT NONE
    
    INTEGER M,P,Q,ierr
    real*8 value,ar(0:PTERM,0:PTERM),ai(0:PTERM,0:PTERM)
    real*8 fv1(0:PTERM),fv2(0:PTERM),fv3(0:PTERM)
    real*8 wr(0:PTERM),wi(0:PTERM),zr(0:PTERM,0:PTERM),zi(0:PTERM,0:PTERM)

!Calculate the eigenvector
    DO 10 M=0,MTERM
        
        IF(DIMAG(BETAK(M)*B)<0.1D0)THEN
            LMAT(0:PTERM,0:PTERM)=CDCOS(BETAK(M)*B)*LMAT1(0:PTERM,0:PTERM)-LMAT2(0:PTERM,0:PTERM)
        ELSE
            LMAT(0:PTERM,0:PTERM)=LMAT1(0:PTERM,0:PTERM)-LMAT2(0:PTERM,0:PTERM)*2.0D0*CDEXP(CI*BETAK(M)*B)/(1.0D0+CDEXP(2.0D0*CI*BETAK(M)*B))
        END IF
        
        
        ar=dreal(LMAT)
        ai=dimag(LMAT)
        CALL cgK(PTERM+1,PTERM+1,ar,ai,wr,wi,1,zr,zi,fv1,fv2,fv3,ierr)
    
        value=cdabs(wr(0)+ci*wi(0))
        do q=1,pterm
            if(cdabs(wr(q)+ci*wi(q))<value)then
                value=cdabs(wr(q)+ci*wi(q))
                p=q
            end if
        enddo
    
        wtp(M,0:PTERM)=zr(0:pterm,p)+ci*zi(0:pterm,p)
    
10  CONTINUE
    
    END
    
    
    