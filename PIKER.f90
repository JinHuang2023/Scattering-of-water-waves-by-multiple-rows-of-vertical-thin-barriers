
    
    SUBROUTINE PIKER_R(N,W)
    
    IMPLICIT NONE
    
    INTEGER N,I,J
    REAL*8 W(n)
    REAL*8 A,C(N)
    
    DO 10 J=2,N
        A=W(J)
        DO 11 I=J-1,1,-1
            IF(W(I)<=A)GOTO 12
            W(I+1)=W(I)
11      CONTINUE
        I=0
12      W(I+1)=A
10  CONTINUE

    END
    
    
    
    SUBROUTINE PIKER(N,W)
    
    IMPLICIT NONE
    
    INTEGER N,I,J
    complex*16 w(n)
    REAL*8 WR(n),WI(N)
    REAL*8 A,B,C(N),D(N)
    
    wr=dimag(w)
    wi=dreal(w)
    
    DO 10 J=2,N
        A=WR(J)
        B=WI(J)
        DO 11 I=J-1,1,-1
            IF(WR(I)<=A)GOTO 12
            WR(I+1)=WR(I)
            WI(I+1)=WI(I)
11      CONTINUE
        I=0
12      WR(I+1)=A
        WI(I+1)=B
10  CONTINUE
    
    w(:)=dcmplx(wi(:),wr(:))
    
    END