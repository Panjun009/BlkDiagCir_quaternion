function [xend,nomres_end,ite,cputime]=fastQuaPCG(Tcol,Trow,b,Ccol,tol,maxiter)

%% Given T which is a Toeplitz quaternion matrix (Tcol and Trow are the 1st column and row respectively), b is a given quaternion vector. This function computes the coefficient vector x, such that x=argmin||b-Ax||
%%  P is preconditioner. For Toeplitz system, we construct circulant matrix as its preconditioner. Here we can use Tony CHan's. 

%% Reference: J. Pan and M. Ng, "Block Diagonalization of Quaternion Circulant Matrices with Applications "   https://arxiv.org/abs/2302.04086  


%%---------- quaternion CG ----------
 
 
m=length(b);
xx0=zerosq(m,1);
 
x0=xx0;
b0=fastQmuliply_Topelvec(Tcol,Trow,x0);
res0=b-b0;  

invcol=fastQinverse_circ(Ccol);
pres0=fastQmuliply_circvec(invcol,res0);
 
dire=pres0;  nores=1e3;

k=0; 
while nores>tol && k<=maxiter
    tic; 
    %dire0=dire;
    bt=fastQmuliply_Topelvec(Tcol,Trow,dire);  % T*dire
    alpha=(res0'*pres0)/(dire'*bt);
    
    x1=x0+alpha*dire;
    
    res1 = res0 - alpha*bt;
    pres1=fastQmuliply_circvec(invcol,res1);
    
    
    nores=norm(res1);
    beta=(res1'*pres1 )/(res0'*pres0);
    dire=pres1+beta*dire;
    
%     dire1=dire;
%     innp=dire0'*A*dire1 
    
    x0=x1;
    res0=res1;
    pres0=pres1;
    
    k=k+1;    
    tim(k)=toc;
end
xend=x1;
nomres_end=nores;
ite=k;
 
cputime=sum(tim);

