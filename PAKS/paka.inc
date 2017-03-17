C
CE    NTOTAL - MAXIMUM TOTAL STORAGE AVAILABLE IN BLANK COMMON A
CS    NTOTAL - MAKSIMALAN PROSTOR U VEKTORU A
C
      integer*8 ntotal
      PARAMETER (ntotal = 400000000)
c academic version
c      PARAMETER (NTOT =  5000000)
c for students
c      PARAMETER (NTOT =  500000)

      REAL A
      COMMON /ZILEGL/ A(ntotal)
c indikator za preko 99999 cvorova
      integer ind56
      common /indikator56/ ind56
        
