////////////////////////////////////////////////////
//MAP431 projet: Estimation d'erreurs a posteriori//
//		Xing WEI & Chongmo LIU		  //
//		     mai,2014			  //
////////////////////////////////////////////////////

//		     Question11.  		  //

clear all
stacksize('max')
//définir des fonction pour résoudre u et simga
function U=solveUh(N,h,f,k,alpha)//pour résoudre uh
    //dimension de U=N
    a=ones(N,1);
    b=ones(N-1,1);
    K=diag(a)*2+diag(b,1)*(-1)+diag(b,-1)*(-1);
    K=sparse(K);
    M=1/6*(diag(a)*4+diag(b,1)+diag(b,-1));
    M=sparse(M);
    Ah=k/h*K+alpha*h*M;
    b=ones(N,1)*h*f;
    U=lusolve(Ah,b);//resoudre le system Ah*u=b
endfunction

function sigma=solveSigmah(N,h,f,k,alpha)//pour résoudre sigma_h
    N=N+2//dimension de sigma=N+2
    a=ones(N,1);
    b=ones(N-1,1);
    K=diag(a)*2+diag(b,1)*(-1)+diag(b,-1)*(-1);
    K(1,1)=1
    K(N,N)=1
    Kp=sparse(K);
    M=(diag(a)*4+diag(b,1)+diag(b,-1));
    M(1,1)=2
    M(N,N)=2
    Mp=sparse(1/6*M)
    Bh=h/k*Mp+1/alpha/h*Kp
    c=zeros(N,1)
    c(1)=1/alpha*f
    c(N)=-1/alpha*f
    sigma=lusolve(Bh,c);//resoudre le system Bh*sigma=c
endfunction

//les paramètres
N=100
h=1/(N+1)
f=1
k=1
alpha=1

Iv=speye(N,N)
Iv=[zeros(1,N);Iv;zeros(1,N)]
Iv=sparse(Iv)//on obtient la matrice creuse Iv

c=[1;1]
//c=sparse(c)
Iw=eye(N+2,N+2)
Iw=kron(Iw,c)
Iw=Iw(2:2*N+3,:)
Iw=sparse(Iw)//on obtient la matrice creuse Iw

D=[-1 1;-1 1]
Dh=speye(N+1,N+1)
Dh=kron(Dh,D)/h//générer la matrice Dh avec la fonction kron
Dh=sparse(Dh)//on obtient la matrice creuse Dh

NN=[2 1;1 2]
NN=sparse(NN)
Nh=speye(N+1,N+1)
Nh=kron(Nh,NN)*h/6
Nh=sparse(Nh)//on obtient la matrice creuse Nh

Fh=f*ones(2*N+2)
Fh=sparse(Fh)//on obtient la matrice creuse Fh

//calculer sigma et u
Uh=solveUh(N,h,f,k,alpha)
Sigmah=solveSigmah(N,h,f,k,alpha)

aa=Iw*Sigmah-k*Dh*Iw*Iv*Uh
bb=Fh-alpha*Iw*Iv*Uh+Dh*Iw*Sigmah
err= 1/k*(Nh*aa)'*aa +1/alpha*(Nh*bb)'*bb //estimation d'erreurs selon l'équation (5)
disp('erreur =')
disp(err)
