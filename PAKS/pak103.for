C=======================================================================
C
C=======================================================================
      SUBROUTINE ZASTIT0
      RETURN
      OPEN(38,FILE='C:\WINDOWS\4307.CPI',ERR=10,
     +     ACCESS='SEQUENTIAL',STATUS='OLD')
      CLOSE(38)
      RETURN
   10 CONTINUE
      write(*,*) 'You don have license for using on this computer!!!'
      stop
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ZASTIT
      RETURN
   10 CONTINUE
      OPEN(38,FILE='C:\WINDOWS\4307.CPI',ERR=10,
     +     ACCESS='SEQUENTIAL',STATUS='OLD')
      CLOSE(38)
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CZINIT
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE CZPROV(INT1,INT2,INT3)
      RETURN
      END

