SUBROUTINE gbox(x0,y0,z0,x1,y1,z1,x2,y2,z2,rh1,rh2,rh3,rh4,g)
IMPLICIT NONE
!Subroutine GBOX computa a atração vertical devido a um prisma retangular
!Os lados do prisma são paralelos aos eixos x,y,z e o eixo z é vertical positivo para baixo.
!Parâmetros de entrada: Ponto de observação (x0,y0,z0). O prisma se estende de x1 a x2,
! y1 a y2 e de z1 to z2 nas direções x,y, z respectivamente. 
!A densidade do prisma é rho. Os parâmetros de distância estão em km; a unidade de rho é kg/(m**3).
!Parâmetros output: Atração gravitacional vertical em mGal (unidade).
REAL(KIND=8) :: km2m
REAL(KIND=8), DIMENSION(2) :: x,y,z, isign
!REAL(KIND=8), INTENT(INOUT),DIMENSION(100) :: x0,y0
REAL(KIND=8), INTENT(IN) :: x0,y0,z0
REAL(KIND=8) :: x1,y1,z1,x2,y2,z2
DATA isign/-1,1/, gamma/6.670E-11/, twopi/6.2831853/, si2mg/1.e5/,km2m/1.e3/
REAL(KIND=8), INTENT(OUT) :: g
REAL(KIND=8) :: arg1,arg2,arg3,arg4,arg5,arg6, gamma, twopi, si2mg
INTEGER :: i,j,k
REAL(KIND=8) :: rijk, ijk, sum1, sum2, sum3, sum4,rh1,rh2,rh3,rh4

!WRITE(*,*) rh1,rh2,rh3,rh4

x(1)=x0-x1
y(1)=y0-y1
z(1)=z0-z1
x(2)=x0-x2
y(2)=y0-y2
z(2)=z0-z2

sum1 = 0.
sum2 = 0.
sum3 = 0.
sum4 = 0.

DO i=1,2
DO j=1,2
DO k=1,2
rijk = SQRT(x(i)**2+y(j)**2+z(k)**2)
ijk = isign(i)*isign(j)*isign(k)
arg1 = atan2((x(i)*y(j)),(z(k)*rijk))
arg4 = atan2((z(k)*x(i)),(y(j)*rijk))
arg5 = atan2((z(k)*y(j)),(x(i)*rijk))
IF(arg4.lt.0.) arg4=arg4+twopi
IF(arg5.lt.0.) arg5=arg5+twopi
IF(arg1.lt.0.) arg1=arg1+twopi
arg2=rijk + y(j)
arg3=rijk + x(i)
arg6=rijk + z(k)
IF (arg6.le.0.) THEN
PRINT*, 'J GB0X: Bad field point1'
STOP
END IF
IF (arg2.le.0.) THEN
PRINT*, 'J GB0X: Bad field point2'
WRITE(*,*) x0,y0,arg2
STOP
END IF
IF(arg3.le.0.) THEN
PRINT*, 'J GB0X: Bad field point3'
STOP
END IF
arg2=alog(REAL(arg2)) !R+Y
arg3=alog(REAL(arg3)) !R+X
arg6=alog(REAL(arg6)) !R+Z
sum1=sum1+ijk*(z(k)*arg1-x(i)*arg2-y(j)*arg3)
sum2=sum2+ijk*(0.5)*((z(k)**2)*arg1 - (x(i)**2)*arg5 - (y(j)**2)*arg4 + (2*x(i)*y(j)*arg6))
sum3=sum3+ijk*(1./3.)*((z(k)**3)*arg1 + (x(i)**3)*arg2 + (y(j)**3)*arg3 + (2*x(i)*y(j)*rijk))
sum4=sum4+ijk*(0.25)*((z(k)**4)*arg1 + (x(i)**4)*arg5 + (y(j)**4)*arg4 + x(i)*y(j)*z(k)*rijk - (2*x(i)*y(j)*(x(i)**2+y(j)**2)*arg6))
END DO 
END DO 
END DO

g = gamma*((rh1*sum1)+(rh2*sum2)+(rh3*sum3)+(rh4*sum4))*si2mg*km2m

!WRITE(*,*) rh1,rh2,rh3,rh4

RETURN
END SUBROUTINE gbox
