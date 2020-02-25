PROGRAM gmodel
IMPLICIT NONE
INTEGER, PARAMETER :: POUT=20, PIN=21, n=10, m=10, GLOCAL=28
INTEGER :: i, j, li, lj
REAL(KIND=8), DIMENSION(0:n) :: x0
REAL(KIND=8), DIMENSION(0:m) :: y0 
REAL(KIND=8) :: xm(0:n), ym(0:m)
REAL(KIND=8) :: z0
REAL(KIND=8) :: x1,y1,z1,x2,y2,z2,rh1,rh2,rh3,rh4,g,g_n,p,q,r,s
REAL(KIND=8) :: xmax,xmin,ymax,ymin,istep,jstep,istepgrav,jstepgrav
REAL(KIND=8), EXTERNAL :: eq
EXTERNAL :: gbox
   

OPEN(UNIT=POUT,FILE='prism.out',STATUS='UNKNOWN') 
OPEN(UNIT=PIN,FILE='para.in',STATUS='UNKNOWN')
OPEN(UNIT=GLOCAL,FILE='glocal.in',STATUS='UNKNOWN')

READ(PIN,*) z0,z1
READ(PIN,*) p,q,r,s
READ(PIN,*) xmin,xmax,ymin,ymax

! Passo para escrever os vetores x() e y()!

istep = (xmax-xmin)/n
jstep = (ymax-ymin)/m
istepgrav = istep/2
jstepgrav = jstep/2

! Gerar a Malha (nxm) - vetores x(n),y(m) !
xm(0) = xmin
ym(0) = ymin
x0(0) = xmin + istepgrav
y0(0) = ymin + jstepgrav

DO i=1,n-1
x0(i)=x0(i-1) + istep
END DO
DO i=1,n 
xm(i)=xm(i-1) + istep
END DO

DO j=1,m-1
y0(j)=y0(j-1) + jstep
END DO
DO j=1,m
ym(j)=ym(j-1) + jstep
END DO

rh1 = p + q*z0 + r*(z0**2) + s*(z0**3)
rh2 = q + 2.*r*z0 + 3.*s*(z0**2)
rh3 = r + 3.*s*z0
rh4 = s

loopcentroy: DO j=0,m-1,1 ! m+1 vetores = m+1 ponto no grid (eixo y)
loopcentrox: DO i=0,n-1,1 ! n+1 vetore = n+1 pontos no grid (eixo x)

g_n = 0

loopyprisma: DO lj=0,m-1,1
loopxprisma: DO li=0,n-1,1

g = 0  ! Entrada de g=0 / cálculo de g para cada prisma
x1 = xm(li)
x2 = xm(li+1)   ! Grafos
y1 = ym(lj)
y2 = ym(lj+1)
z2 = eq(x0(li),y0(lj)) ! z2 centrado no prisma utilizando equação modelo
CALL gbox(x0(i),y0(j),z0,x1,y1,z1,x2,y2,z2,rh1,rh2,rh3,rh4,g) ! Subroutina para cálculo de g de cada prisma
g_n=g_n+g ! soma da gravidade em um ponto devido a contribuição de todos os prisma


END DO loopxprisma
END DO loopyprisma

WRITE(POUT,*) x0(i), y0(j), -eq(x0(i),y0(j))
WRITE(GLOCAL,*) x0(i), y0(j), g_n

END DO loopcentrox  
END DO loopcentroy


CLOSE(POUT)
CLOSE(PIN)
CLOSE(GLOCAL)

END PROGRAM gmodel

INCLUDE 'gbox.f90'

REAL(KIND=8) FUNCTION eq(x,y)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: x,y

!eq=(sin((y/2)*cos(x/2)) -  sin(x)*cos(y)) + 2
eq= 0.012*(x**2) + 0.0012*(y**2)
!eq= -0.012*(x**2 + y**2) + 0.12*(x + y)  

END FUNCTION eq
