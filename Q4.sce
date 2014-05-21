////////////////////////////////////////////////////
//MAP431 projet: Estimation d'erreurs a posteriori//
//		Xing WEI & Chongmo LIU		  //
//		     mai,2014			  //
////////////////////////////////////////////////////

//		     Question4.  		  //

clear all
//les paramètres
N=100;
h=1/(N+1)
f=1;
k=1;
alpha=1;

a=ones(N,1);
b=ones(N-1,1);
K=diag(a)*2+diag(b,1)*(-1)+diag(b,-1)*(-1);//on génère la matrice K
K=sparse(K);//on convertie la matrice K en matrice creuse
M=1/6*(diag(a)*4+diag(b,1)+diag(b,-1));
M=sparse(M);//on convertie la matrice M en matrice creuse
Ah=k/h*K+alpha*h*M;
b=ones(N,1)*h*f;

u=lusolve(Ah,b);//resoudre le system Ah*u=b

xtitle('graph de u','x','u(x)')
plot(linspace(0,1,N),u)//tracer lacourbe de u
