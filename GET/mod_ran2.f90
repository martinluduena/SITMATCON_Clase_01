MODULE mod_ran2
contains
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL(8) ran2,AMM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AMM=1./IM1,IMM1=IM1-1,&
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/24453211/, iv/NTAB*0/, iy/0/
      IF(idum.le.0)THEN
      idum=max(-idum,1)
      idum2=idum
      DO 11 j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      IF(idum.lt.0) idum=idum+IM1
      IF(j.le.NTAB) iv(j)=idum
11    CONTINUE
      iy=iv(1)
      END IF
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      IF(idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      IF(idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      IF(iy.lt.1)iy=iy+IMM1
      ran2=min(AMM*iy,RNMX)
      RETURN
      END FUNCTION ran2

ENDMODULE mod_ran2

