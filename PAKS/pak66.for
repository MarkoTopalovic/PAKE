C=======================================================================
C
C=======================================================================
      SUBROUTINE SNEL6(S,FE,DUZ,POVR,ZI,YI,TK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   NELINEARNA MATRICA KRUTOSTI TANKOZIDNE GREDE  SN
C
      DIMENSION S(12,*),FE(*)
      DUZ2=DUZ*DUZ
      DUZ3=DUZ2*DUZ
      DUZ30=30.D0*DUZ
      PN=0.5D0*(FE(7)-FE(1))
      A=PN/DUZ
      B=12.D0*PN*ZI/POVR/DUZ3+36.D0*PN/DUZ30
      C=12.D0*PN*YI/POVR/DUZ3+36.D0*PN/DUZ30
      D=PN*TK/POVR/DUZ2
      E=4.D0*PN*YI/POVR/DUZ+2.D0*PN*DUZ/15.D0
      F=4.D0*PN*ZI/POVR/DUZ+2.D0*PN*DUZ/15.D0
      G=FE(5)/DUZ
      H=FE(6)/DUZ
      G1=FE(11)/DUZ
      H1=FE(12)/DUZ
      AK=6.D0*PN*ZI/POVR/DUZ2+0.1D0*PN
      AL=6.D0*PN*YI/POVR/DUZ2+0.1D0*PN
      AM=(FE(6)+FE(12))/6.D0
      AN=(FE(5)+FE(11))/6.D0
      E1=2.D0*PN*YI/POVR/DUZ-PN*DUZ/30.D0
      F1=2.D0*PN*ZI/POVR/DUZ-PN*DUZ/30.D0
C...  RASPOREDJIVANJE CLANOVA U MATRICU S = KL + KNL
      S(1,1)=S(1,1)+A      
      S(1,5)=S(1,5)+G      
      S(1,6)=S(1,6)+H      
      S(2,2)=S(2,2)+B      
      S(2,4)=S(2,4)-G      
      S(2,6)=S(2,6)+AK      
      S(3,3)=S(3,3)+C      
      S(3,4)=S(3,4)-H      
      S(3,5)=S(3,5)-AL      
      S(4,4)=S(4,4)+D      
      S(4,5)=S(4,5)-AM      
      S(4,6)=S(4,6)+AN      
      S(5,5)=S(5,5)+E      
      S(6,6)=S(6,6)+F
C
      S(7,7)=S(7,7)+A      
      S(7,11)=S(7,11)-G1      
      S(7,12)=S(7,12)-H1      
      S(8,8)=S(8,8)+B      
      S(8,10)=S(8,10)+G1      
      S(8,12)=S(8,12)-AK      
      S(9,9)=S(9,9)+C      
      S(9,10)=S(9,10)+H1      
      S(9,11)=S(9,11)+AL      
      S(10,10)=S(10,10)+D      
      S(10,11)=S(10,11)+AM      
      S(10,12)=S(10,12)-AN      
      S(11,11)=S(11,11)+E      
      S(12,12)=S(12,12)+F
C      
      S(1,7)=S(1,7)-A      
      S(1,11)=S(1,11)+G1      
      S(1,12)=S(1,12)+H1      
      S(2,8)=S(2,8)-B      
      S(2,10)=S(2,10)-G1      
      S(2,12)=S(2,12)+AK      
      S(3,9)=S(3,9)-C      
      S(3,10)=S(3,10)-H1      
      S(3,11)=S(3,11)-AL      
      S(4,8)=S(4,8)+G      
      S(4,9)=S(4,9)+H      
      S(4,10)=S(4,10)-D      
      S(4,11)=S(4,11)-AM      
      S(4,12)=S(4,12)+AN      
      S(5,7)=S(5,7)-G      
      S(5,9)=S(5,9)+AL      
      S(5,10)=S(5,10)+AM      
      S(5,11)=S(5,11)+E1      
      S(6,7)=S(6,7)-H
      S(6,8)=S(6,8)-AK
      S(6,10)=S(6,10)-AN
      S(6,12)=S(6,12)+F1
C
      DO 10 I=1,12
      DO 10 J=I,12
   10 S(J,I)=S(I,J)
      RETURN
      END
