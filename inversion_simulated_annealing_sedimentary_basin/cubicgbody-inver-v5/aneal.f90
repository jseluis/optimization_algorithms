PROGRAM aneal
IMPLICIT NONE
! tfactor = fator de resfriamento / redução da temperatura
! nsteps = numero de passos em um bloco
INTEGER, PARAMETER :: ndim=1
INTEGER, PARAMETER :: ENOUT=11, XOUT=12
INTEGER :: nsteps=100000000, istepsim, idum
INTEGER :: naccepted
REAL(KIND=16) :: fx, temperatura, tfactor, t1, t2
REAL(KIND=16), EXTERNAL ::  func, ran3
REAL(KIND=16), DIMENSION(ndim)  :: x, xsum, x2sum, xmin, y, ymin
REAL(KIND=16)  :: ensum, en2sum, enav, cv !computar o calor específico
REAL(KIND=16)  :: aratio, anorm, emin, sigmax, z, t
EXTERNAL step

OPEN(UNIT=ENOUT,FILE="en.out",STATUS="UNKNOWN")
OPEN(UNIT=XOUT,FILE="x.out",STATUS="UNKNOWN")

WRITE(*,*) 'Chute um valor inicial para x,y,z,t:'
READ(*,*) x,y,z,t
fx=func(x,y,z,t)
emin=fx
naccepted=0

WRITE(11,'("iniciando com x=",5f20.10)') x
WRITE(12,'("valor da função custo f=",5f20.10)') fx

WRITE(*,*) 'Entrada inicial de (alta) temperatura (e.g., 10)'
10 READ(*,*) temperatura
  IF(temperatura.lt.0.001) THEN
     temperatura=0.001
     WRITE(*,*) '  resete a temperatura =', temperatura
  ELSEIF (temperatura.gt.10) THEN
  WRITE(*,*) temperatura,'é uma Temperatura muito alta'
  GOTO 10
  ENDIF


!tfactor = 0.9d0     Fator de cozimento

WRITE(*,*) '  Entrada do fator de cozimento(e.g., 0.9)'
READ(*,*) tfactor
WRITE(*,*) '  Número de passos por bloco (equilibrio, e.g., 1000)'
READ(*,*) nsteps

! Gerador aleatorio

WRITE(*,*) '  Leitura do idum = gerador aleatorio, e.g., 18'
READ(*,*) idum

!Antes de iniciar o Loop, termalizar o bloco
DO istepsim=1,nsteps
  CALL step(x(1:1),fx,temperatura,func,naccepted,idum,y,z,t,nsteps)
END DO

WRITE(*,'("#Termalizado para",i5," steps.")') naccepted
WRITE(*,'("#Loop com",i5," steps.")') nsteps


WRITE(ENOUT,'(a8,2x,a11,3(2x,a8),a11)') &
       "#temp ","<E>  ", "C_v  ", "sigma_x  ", "ac.ratio","E_min"
  
WRITE(XOUT,'("#temp, ndim colunas para <x>, ndim colunas para x_min, <y>, y_min")')

CALL cpu_time(t1)
  
ciclo_anneal: DO
    
    !Primeiro bloco de um passo metropolis
    naccepted = 0 
    !Acompanhamos <E> e <E^2> para calcular C_v
    ensum=0d0; en2sum=0d0; xsum=0d0; x2sum=0d0
    DO istepsim=1,nsteps
      CALL step(x(1:ndim),fx,temperatura,func,naccepted,idum,y,z,t,nsteps)
      ensum = ensum+fx
      en2sum = en2sum+fx**2
      xsum = xsum + x
      x2sum = x2sum + x**2
    IF (fx<emin) THEN
        emin = fx
        xmin = x
        ymin = y
    END IF
    END DO

   anorm=1d0/FLOAT(nsteps)
    enav = ensum*anorm
    cv = (en2sum*anorm-enav**2)/temperatura**2
    sigmax = SQRT(SUM( x2sum*anorm - (xsum*anorm)**2 ))
    aratio=naccepted*anorm

    WRITE(ENOUT,'(e8.2,2x,e11.5,3(2x,e8.2),2x,e11.5)'), temperatura, enav, cv, sigmax, aratio, emin 

    WRITE(XOUT,'(20(e22.15,2x))'), temperatura, x(1:ndim), xmin(1:ndim),y(1:ndim), ymin(1:ndim)

    !Resfriamento do sistema (decaimento da temperatura)
    temperatura = temperatura * tfactor
   
    IF ((emin < 10d-8).and.(temperatura < 10d-10)) EXIT ciclo_anneal
  END DO ciclo_anneal

CALL cpu_time( t2 )
WRITE(*,*) 'FINISH-start', t2-t1

CLOSE(XOUT)
CLOSE(ENOUT)
END PROGRAM aneal

!RAN3 - Função gerador aleatório

REAL(KIND=16) FUNCTION ran3(idum)
!Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value
!to initialize or reinitialize the sequence.
INTEGER idum
INTEGER MBIG,MSEED,MZ
!REAL MBIG,MSEED,MZ
REAL(KIND=16) FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
!According to Knuth, any large mbig, and any smaller (but still large) mseed can be substituted
!for the above values.

INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55) !The value 55 is special and should not be modified; see
!REAL mj,mk,ma(55)            !Knuth.
SAVE iff,inext,inextp,ma
DATA iff /0/
    
    IF(idum.lt.0.or.iff.eq.0) THEN !Initialization.
iff=1
mj=abs(MSEED-abs(idum)) !Initialize ma(55) using the seed idum and the large number mseed
mj=mod(mj,MBIG)                                
ma(55)=mj
mk=1
DO i=1,54     !Now initialize the rest of the table,
ii=mod(21*i,55) !in a slightly random order,
ma(ii)=mk !with numbers that are not especially random.
mk=mj-mk
IF(mk.lt.MZ)mk=mk+MBIG
END DO

DO  k=1,4 !We randomize them by “warming up the generator.”
DO  i=1,55
ma(i)=ma(i)-ma(1+mod(i+30,55))
IF(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
ENDDO 
ENDDO 
inext=0 !Prepare indices for our first generated number.
inextp=31 !The constant 31 is special; see Knuth.
idum=1
ENDIF
inext=inext+1 !Here is where we start, except on initialization. Increment
IF(inext.eq.56)inext=1 !inext, wrapping around 56 to 1.
inextp=inextp+1 !Ditto for inextp.
IF(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp) !Now generate a new random number subtractively.
IF(mj.lt.MZ)mj=mj+MBIG !Be sure that it is in range.
ma(inext)=mj !Store it,and output the derived uniform deviate.
ran3=mj*FAC 
return
END

INCLUDE 'step.f90'	

! Função Minimizada (custo)
REAL(KIND=16) FUNCTION func(x,y,z,t)
  IMPLICIT NONE
  REAL(KIND=16), INTENT(IN), DIMENSION(1) :: x,y
  REAL(KIND=16) :: z,t
  func = SQRT((x(1)*t + (y(1)*(t**2))/2 - z)**2)
END FUNCTION



