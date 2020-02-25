SUBROUTINE step(xm,ym,x0,y0,z0,z1,z,fx,temperatura,func,naccepted,idum,g_campo,nsteps,p1,q1,r1,s1,g_calc,emin)
IMPLICIT NONE
  !Executar o passo metropolis
  !Variáveis e argumentos
  REAL(KIND=8), INTENT(INOUT), DIMENSION(0:9,0:9) :: g_campo,g_calc
  REAL(KIND=8), INTENT(INOUT), DIMENSION(0:9,0:9) :: z
  REAL(KIND=8), INTENT(INOUT) :: fx,p1,q1,r1,s1
  REAL(KIND=8), INTENT(IN) :: temperatura,z0,z1,emin
  INTEGER, INTENT(INOUT) :: naccepted
  INTEGER, INTENT(IN) :: idum
!  REAL(KIND=8), INTENT(IN), DIMENSION(:) :: x0
!  REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xm ! n-1 posições de centros dos prismas (eixo x)
!  REAL(KIND=8), INTENT(IN), DIMENSION(:) :: y0
!  REAL(KIND=8), INTENT(IN), DIMENSION(:) :: ym ! m-1 posições de centros dos prismas (eixo y) 
  REAL(KIND=8), INTENT(IN) :: x0(0:9),xm(0:10),y0(0:9),ym(0:10)
  !Variáveis locais
  REAL(KIND=8) :: fx_n,g_n,g,x1,x2,y1,y2,z2,c,rh1,rh2,rh3,rh4,p_n,q_n,r_n,s_n
  REAL(KIND=8), DIMENSION(0:9,0:9) :: z_n, g_sub
  REAL(KIND=8) :: e=0.2   ! Escala escolhida para FUNC específica
  INTEGER :: a,nsteps,i,j,li,lj
  INTEGER, PARAMETER :: nn=10, mm=10, qq=9, ww=9
  REAL(KIND=8), EXTERNAL :: ran3, func
  EXTERNAL gbox

IF ((emin<2).and.(emin>1)) THEN
e=0.1
ELSE IF ((emin<1).and.(emin>0.1)) THEN
e=0.05
ELSE IF((emin<0.1).and.(emin>0.01)) THEN
e=0.01
ELSE
END IF

!Gerar novas alturas aleatoriamente
DO j=0,ww
DO i=0,qq
z_n(i,j) = REAL(z(i,j) + 2*e*z(i,j)*(ran3(idum) - 0.5)) ! movimento estocático no alcance
!  write(*,*) z_n(i,j),z(i,j)
END DO 
END DO
!FUNCIONANDO - DEBUGADO

!Gerar novos parâmetros do modelo aleatoriamente
p_n = REAL(p1 + 2*e*p1*(ran3(idum) - 0.5))
q_n = REAL(q1 + 2*e*q1*(ran3(idum) - 0.5))
r_n = REAL(r1 + 2*e*r1*(ran3(idum) - 0.5))  
s_n = REAL(s1 + 2*e*s1*(ran3(idum) - 0.5)) 
!Gerar densidade com base nos novos parâmetros do modelo 
rh1 = p_n + q_n*z0 + r_n*(z0**2) + s_n*(z0**3)
rh2 = q_n + 2.*r_n*z0 + 3.*s_n*(z0**2)
rh3 = r_n + 3.*s_n*z0
rh4 = s_n

!Calcular a anomalia devido às novas alturas
  
loopcentroy: DO j=0,ww,1 ! m+1 vetores = m+1 ponto no grid (eixo y)
loopcentrox: DO i=0,qq,1 ! n+1 vetore = n+1 pontos no grid (eixo x)

g_n = 0

loopyprisma: DO lj=0,ww,1
loopxprisma: DO li=0,qq,1

g = 0  ! Entrada de g=0 / cálculo de g para cada prisma
x1 = xm(li)
x2 = xm(li+1)   ! Grafos
y1 = ym(lj)
y2 = ym(lj+1)
z2 = z_n(li,lj)
CALL gbox(x0(i),y0(j),z0,x1,y1,z1,x2,y2,z2,rh1,rh2,rh3,rh4,g) ! Subroutina para cálculo de g de cada prisma
g_n=g_n+g ! soma da gravidade em um ponto devido a contribuição de todos os prisma

END DO loopxprisma
END DO loopyprisma

g_sub(i,j)= g_n

!WRITE(*,*) g_campo(i,j),g_calc(i,j),g_sub(i,j)
!FUNCIONANDO - DEBUGADO
END DO loopcentrox  
END DO loopcentroy

! novo valor da função custo/objetiva

fx_n=0    
c=0
DO j=0,ww
DO i=0,qq
c = func(g_campo(i,j),g_sub(i,j))
fx_n= c + fx_n
END DO 
END DO

a=nsteps/3 ! número de aceitação de passos ruins!

 externo: IF (fx_n < fx) THEN
  
  fx = fx_n                ! Condicional salvo com sucesso
  DO j=0,ww
  DO i=0,qq      
  
 z(i,j)= z_n(i,j)
 g_calc(i,j)=g_sub(i,j) ! Atualizar todas as alturas para o programa principal
  
  END DO 
  END DO
  
p1=p_n
q1=q_n
r1=r_n
s1=s_n

 ELSE 
  interno : IF ((EXP(-(fx_n - fx)/temperatura) > ran3(idum)) .and. (naccepted < a)  ) THEN
   
  fx = fx_n
  p1=p_n
  q1=q_n
  r1=r_n
  s1=s_n                ! Condicional salvo com sucesso
  
  DO j=0,ww
  DO i=0,qq      
  z(i,j)= z_n(i,j) ! Atualizar todas as alturas para o programa principal
  END DO 
  END DO
  
  naccepted = naccepted+1
  END IF interno
 
  END IF externo

END SUBROUTINE step
