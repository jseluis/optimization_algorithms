PROGRAM funcao
IMPLICIT NONE
INTEGER, PARAMETER :: gsaida=23, GLOCAL=24, q=9, w=9
INTEGER :: i,j,n
REAL(KIND=16), DIMENSION(0:q,0:w) :: g_campo, g_calc ! 
REAL(KIND=16) :: custo ,fx
REAL(KIND=16), EXTERNAL :: func

OPEN(UNIT=gsaida,FILE='gsaida.out',STATUS='UNKNOWN')
OPEN(UNIT=GLOCAL,FILE='glocal.in',STATUS='UNKNOWN')

n=0
custo = 0 
DO j=0,w
DO i=0,q
READ(GLOCAL,*) g_campo(i,j)
READ(gsaida,*) g_calc(i,j)
fx = func(g_campo(i,j),g_calc(i,j))
custo=custo + fx
n=n+1
END DO 
END DO

WRITE(*,*) 'Saida',n, custo

CLOSE(gsaida)
CLOSE(GLOCAL)
STOP
END PROGRAM funcao


REAL(KIND=16) FUNCTION func(g,gc) 
IMPLICIT NONE
REAL(KIND=16), INTENT(IN) :: g,gc
    func = SQRT((g - gc)**2)

END FUNCTION func
