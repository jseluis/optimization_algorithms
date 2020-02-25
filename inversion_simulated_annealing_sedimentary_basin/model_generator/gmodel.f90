PROGRAM gmodel
IMPLICIT NONE
INTEGER, PARAMETER :: POUT=20, PIN=21,  n=100, m=100, GLOCAL=28
INTEGER :: i, j 
REAL(KIND=8), DIMENSION(0:n) :: x0
REAL(KIND=8), DIMENSION(0:m) :: y0 
REAL(KIND=8) :: z0
REAL(KIND=8) :: x1,y1,z1,x2,y2,z2,rh1,rh2,rh3,rh4,g,p,q,r,s
REAL(KIND=8) :: xmax,xmin,ymax,ymin,istep,jstep,istepgrav,jstepgrav
EXTERNAL :: gbox
   

OPEN(UNIT=POUT,FILE='prism.out',STATUS='UNKNOWN') 
OPEN(UNIT=PIN,FILE='para.in',STATUS='UNKNOWN')
OPEN(UNIT=GLOCAL,FILE='glocal.in',STATUS='UNKNOWN')

READ(PIN,*) x1,y1,z1,x2,y2,z2
READ(PIN,*) z0,z1
READ(PIN,*) p,q,r,s
READ(PIN,*) xmin,xmax,ymin,ymax

!WRITE(*,*) 'Escreva: x1,y1,z1,x2,y2,z2,p,q,r,s)'
!READ(*,*) x1, y1, z1, x2, y2, z2, p, q, r, s
!WRITE(*,*) 'Escreva: xmin,xmax,ymin,ymax'
!READ(*,*) xmin, xmax, ymin, ymax


! Passo para escrever os vetores x() e y()!

istep = (xmax-xmin)/n
jstep = (ymax-ymin)/m
istepgrav = istep/2 
jstepgrav = jstep/2

! Gerar a Malha (nxm) - vetores x(n),y(m) !

x0(0) = xmin + istepgrav
y0(0) = ymin + jstepgrav

DO i=1,n-1
x0(i)=x0(i-1) + istep
!WRITE(POUT,*) xi(i)
END DO

DO j=1,m-1
y0(j)=y0(j-1) + jstep
!WRITE(POUT,*) y0(j)
END DO

! Posição do prisma, densidade
!WRITE(PIN,*) 'x1','y1','z1','x2','y2','z2'
!WRITE(PIN,*) x1,y1,z1,x2,y2,z2
!WRITE(PIN,*) 'rho','istep','jstep','xmin','xmax','ymin','ymax'
!WRITE(PIN,*) istep,jstep,xmax,xmin,ymax,ymin
!WRITE(POUT,'(2x,5a)') 'x0:            ','y0:               ', 'g(mGal):'

rh1 = p + q*z0 + r*(z0**2) + s*(z0**3)
rh2 = q + 2.*r*z0 + 3.*s*(z0**2)
rh3 = r + 3.*s*z0
rh4 = s

DO j=0,m-1,1
DO i=0,n-1,1
CALL gbox(x0(i),y0(j),z0,x1,y1,z1,x2,y2,z2,rh1,rh2,rh3,rh4,g)
WRITE(POUT,*) x0(i), y0(j), g
WRITE(GLOCAL,*) g
END DO  
END DO 

CLOSE(POUT)
CLOSE(PIN)
CLOSE(GLOCAL)

END PROGRAM gmodel

INCLUDE 'gbox.f90'
