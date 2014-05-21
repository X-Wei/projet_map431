////////////////////////////////////////////////////
//MAP431 projet: Estimation d'erreurs a posteriori//
//		Xing WEI & Chongmo LIU		  //
//		     mai,2014			  //
////////////////////////////////////////////////////

//		     Question6.  		  //

clear all
//les paramètres
N=100;
h=1/(N+1)
f=1;
k=1;
alpha=1;

N=N+2//dimension de Wh=N+2
a=ones(N,1);
b=ones(N-1,1);
K=diag(a)*2+diag(b,1)*(-1)+diag(b,-1)*(-1);
K(1,1)=1
K(N,N)=1//on génère la matrice K
Kp=sparse(K)//convertir K en matrice creuse Kp 
M=(diag(a)*4+diag(b,1)+diag(b,-1));
M(1,1)=2
M(N,N)=2
Mp=sparse(1/6*M)//générer la matrice creuse Mp

Bh=h/k*Mp+1/alpha/h*Kp//on obtient la matrice Bh

c=zeros(N,1)
c(1)=1/alpha*f
c(N)=-1/alpha*f

sigma=lusolve(Bh,c);//resoudre le système Bh*sigma=c

xtitle('graph de delta_h','x','delta(x)')
plot(linspace(0,1,N),sigma)//tracer le flux sigma
