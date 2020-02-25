PROGRAM prism
IMPLICIT NONE
INTEGER, PARAMETER :: POUT=20, PIN=21, gsaida=23, PINN=25, GLOCAL=24, n=10, m=10, q=9, w=9
INTEGER :: i,j,li,lj
!REAL(KIND=8), DIMENSION(0:n) :: x0
!REAL(KIND=8), DIMENSION(0:q) :: xm ! n-1 posições de centros dos prismas (eixo x)
!REAL(KIND=8), DIMENSION(0:m) :: y0
!REAL(KIND=8), DIMENSION(0:w) :: ym ! m-1 posições de centros dos prismas (eixo y)
REAL(KIND=8) :: x0(0:q),xm(0:n),y0(0:w),ym(0:m)
REAL(KIND=8), DIMENSION(0:q,0:w) :: z, g_campo, g_calc, z_min ! matriz de alturas e gravidade medida
REAL(KIND=8) :: z0,z1,z2,x1,x2,y1,y2
REAL(KIND=8) :: g,g_n,p1,q1,r1,s1,rh1x,rh2x,rh3x,rh4x !rho = densidade do prisma 
REAL(KIND=8) :: xmax,xmin,ymax,ymin,istep,jstep,istepgrav,jstepgrav
REAL(KIND=8), PARAMETER :: gamma=6.670E-11, twopi=6.2831853
! Inversão Utilizando o Simulated Annealing
! tfactor = fator de resfriamento / redução da temperatura
! nsteps = numero de passos em um bloco
INTEGER, PARAMETER :: ENOUT=11, XOUT=12
INTEGER :: nsteps=100000000, istepsim, idum
INTEGER :: naccepted
REAL(KIND=8) :: fx, temperatura, tfactor, t1, t2, custo
REAL(KIND=8), EXTERNAL ::  func, ran3
REAL(KIND=8)  :: emin,p1min,q1min,r1min,s1min
EXTERNAL step,termal

OPEN(UNIT=POUT,FILE='prism.out',STATUS='UNKNOWN') 
OPEN(UNIT=PIN,FILE='para.out',STATUS='UNKNOWN')
OPEN(UNIT=PINN,FILE='para.in',STATUS='UNKNOWN')
OPEN(UNIT=ENOUT,FILE="en.out",STATUS="UNKNOWN")
OPEN(UNIT=XOUT,FILE="x.out",STATUS="UNKNOWN")
OPEN(UNIT=GLOCAL,FILE="glocal.in",STATUS="UNKNOWN")
OPEN(UNIT=gsaida,FILE="gsaida.out",STATUS="UNKNOWN")

!temperatura = alta temperatura (e.g., 10)
!tfactor = 0.9d0 - Fator de resfriamento
!nsteps = Número de passos por bloco (equilibrio, e.g., 1000)
!idum = Gerador aleatorio , e.g., 18

READ(PINN,*) z0,z1
READ(PINN,*) p1,q1,r1,s1
READ(PINN,*) xmin,xmax,ymin,ymax
READ(PINN,*) temperatura,tfactor,nsteps,idum

WRITE(PIN,*) 'Parâmetros iniciais da inversão:'
WRITE(PIN,*) 'p1 = ',p1
WRITE(PIN,*) 'q1 = ',q1
WRITE(PIN,*) 'r1 = ',r1
WRITE(PIN,*) 's1 = ',s1
WRITE(PIN,*) 'z0 = ',z0
WRITE(PIN,*) 'z1 = ',z1
WRITE(PIN,*) 'xmin = ',xmin,' xmax = ',xmax
WRITE(PIN,*) 'ymin = ',ymin,' ymax = ',ymax
WRITE(PIN,*) 'temperatura =',temperatura
WRITE(PIN,*) 'tfactor =',tfactor
WRITE(PIN,*) 'nsteps =',nsteps
WRITE(PIN,*) 'idum =',idum

! Passo para escrever os vetores x() e y()!
istep = (xmax-xmin)/n
jstep = (ymax-ymin)/m
istepgrav = 0.5*istep ! n=número de blocos eixo x
jstepgrav = 0.5*jstep ! m = número de blocos eixo y e x*y = Ntotal de blocos.

! Gerar a Malha (nxm) - Prismas e pontos centrados na superfície !
xm(0) = xmin
ym(0) = ymin
x0(0) = xmin + istepgrav  
y0(0) = ymin + jstepgrav

!q=n-1 ! número de pontos na superfície centrados nos prismas
!w=m-1 ! número de pontos na superfície centrados nos prismas

! Leitura dos pontos dos prismas (eixo x)
DO i=1,n 
xm(i)=xm(i-1) + istep
END DO
! Leitura dos pontos centrados nos prismas (eixo x)
DO i=1,q
x0(i)=x0(i-1) + istep
END DO
! Leitura dos pontos dos prismas (eixo y)
DO j=1,m
ym(j)=ym(j-1) + jstep
END DO
! Leitura dos pontos centrados nos prismas (eixoy)
DO j=1,w
y0(j)=y0(j-1) + jstep
END DO
WRITE(PIN,*) 'istep = ',istep
WRITE(PIN,*) 'jstep = ',jstep

!Calcular rh1,rh2,rh3,rh4 inicial

rh1x = p1 + q1*z0 + r1*(z0**2) + s1*(z0**3)
rh2x = q1 + 2.*r1*z0 + 3.*s1*(z0**2)
rh3x = r1 + 3.*s1*z0
rh4x = s1

! Leitura do GLOCAL e chute inicial das alturas

DO j=0,w
DO i=0,q
READ(GLOCAL,*) g_campo(i,j)
z(i,j) = g_campo(i,j)/(gamma*twopi*rh1x*1.e8)
!WRITE(gsaida,*) x0(i),y0(j),g_campo(i,j),z(i,j)
END DO 
END DO

loopcentroy: DO j=0,w,1 ! m+1 vetores = m+1 ponto no grid (eixo y)
loopcentrox: DO i=0,q,1 ! n+1 vetore = n+1 pontos no grid (eixo x)

g_n = 0

loopyprisma: DO lj=0,w,1
loopxprisma: DO li=0,q,1

g = 0  ! Entrada de g=0 / cálculo de g para cada prisma
x1 = xm(li)
x2 = xm(li+1)   ! Grafos
y1 = ym(lj)
y2 = ym(lj+1)
z2 = z(li,lj)
CALL gbox(x0(i),y0(j),z0,x1,y1,z1,x2,y2,z2,rh1x,rh2x,rh3x,rh4x,g) ! Subroutina para cálculo de g de cada prisma
g_n=g_n+g ! soma da gravidade em um ponto devido a contribuição de todos os prisma

END DO loopxprisma
END DO loopyprisma

g_calc(i,j)= g_n
!WRITE(gsaida,*) x0(i),y0(j),g_campo(i,j),g_calc(i,j),z(i,j)
END DO loopcentrox  
END DO loopcentroy

!! Termalizar o sistema

fx = 0 
DO j=0,w
DO i=0,q
custo = func(g_campo(i,j),g_calc(i,j)) ! Calculo da função custo inicial
fx = fx + custo
END DO 
END DO

emin=fx ! salva a função custo como emin
naccepted=0

WRITE(*,*) 'Cost Function: ',fx

!Antes de iniciar o Loop, termalizar o bloco
!!DEBUGUEI ATÉ AQUI!

ll: DO istepsim=1,nsteps

CALL termal(xm(0:n),ym(0:m),x0(0:q),y0(0:w),z0,z1,z(0:q,0:w),fx,temperatura,func,naccepted,idum,g_campo(0:q,0:w),nsteps, &
p1,q1,r1,s1,g_calc(0:q,0:w))
WRITE(*,*) fx
IF(fx<8) THEN
EXIT ll
END IF
END DO ll

DO j=0,w
DO i=0,q
!WRITE(gsaida,*) x0(i),y0(j),z(i,j),g_campo(i,j),g_calc(i,j) ! Calculo da função custo inicial
END DO 
END DO

emin=fx ! salva a função custo como emin

WRITE(*,*) 'Saida-termalizado',fx


WRITE(*,'("#Termalizado para",i5," steps.")') naccepted

!! Fim da termalização

! Início do algoritmo - Simulated Annealing

CALL cpu_time(t1)
  
ciclo_anneal: DO
    
    !Primeiro bloco de um passo metropolis
   
naccepted = 0 
    
    
DO istepsim=1,nsteps
    
CALL step(xm(0:n),ym(0:m),x0(0:q),y0(0:w),z0,z1,z(0:q,0:w),fx,temperatura,func,naccepted,idum,g_campo(0:q,0:w),nsteps, &
p1,q1,r1,s1,g_calc(0:q,0:w),emin)
    
IF (fx<emin) THEN
emin = fx
p1min = p1
q1min = q1
r1min = r1
s1min = s1
WRITE(*,*) emin,p1min,q1min,r1min,s1min

DO j=0,w
DO i=0,q
z_min(i,j)=z(i,j) ! grava o z_min (fora do mínimo local)
END DO 
END DO	

END IF
  
END DO

WRITE(ENOUT,*) temperatura, emin 

temperatura = temperatura * tfactor

IF (emin < 5) THEN
DO j=0,w
DO i=0,q
WRITE(XOUT,*) x0(i),y0(j),z_min(i,j)  
END DO 
END DO
END IF

IF (emin < 4) THEN
DO j=0,w
DO i=0,q
WRITE(XOUT,*) x0(i),y0(j),z_min(i,j)  
END DO 
END DO
EXIT ciclo_anneal
END IF

IF ((temperatura < 0.001).and.(emin > 0.1)) THEN
temperatura = 7.0
END IF
!Resfriamento do sistema (decaimento da temperatura)


END DO ciclo_anneal

CALL cpu_time( t2 )
WRITE(PIN,*) 'FINISH-start', t2-t1

CLOSE(XOUT)
CLOSE(ENOUT)
CLOSE(POUT)
CLOSE(PIN)
CLOSE(PINN)
CLOSE(GLOCAL)
CLOSE(gsaida)

END PROGRAM prism 

!RAN3 - GERADOR ALEATORIO
REAL(KIND=8) FUNCTION ran3(idum)
!Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value
!to initialize or reinitialize the sequence.
INTEGER idum
INTEGER MBIG,MSEED,MZ
!REAL(KIND=8) MBIG,MSEED,MZ
REAL(KIND=8) :: FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
!According to Knuth, any large mbig, and any smaller (but still large) mseed can be substituted
!for the above values.

INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55) !The value 55 is special and should not be modified; see
!REAL(KIND=8) mj,mk,ma(55)            !Knuth.
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
RETURN
END

! Função Minimizada (custo)

REAL(KIND=8) FUNCTION func(g,gc)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: g, gc
    func = SQRT((g - gc)**2)
END FUNCTION

INCLUDE 'gbox.f90'

INCLUDE 'step.f90'

INCLUDE 'termal.f90'
