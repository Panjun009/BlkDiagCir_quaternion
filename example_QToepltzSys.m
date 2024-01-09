%% the code shows an example solving the AR(1) and AR(2) system by fast circulant preconditioner method


 clc
clear all
addpath('qtfm')
% 
% %% ======================== generate model by AR (1) ====================

m = 200 ; % sample number
eta=0.3; pho=0.5; x1=randq(1); xL=length(x1);
[X,y0]=QuaAR1(x1,pho,eta,m);

numbl=4;
n=m/numbl; % window width
% 
%% ======================== generate model by AR (2) ====================

 % n=100;% window width
 % numbl=2;
 % 
 %  m = n*numbl; % sample number
 % eta=0.025; toa=0.9; tob=0.1; x1=randq(1); x2=randq(1);xL=length(x1);
 % [X,y0]=QuaAR2(x1,x2,toa,tob,eta,m);
  
 
%% +++++++++++ by correlation method (T1) Note T1 is Toeplitz +++++++++++

d=(n-1)*xL;
T1 = windowMat(X,n,1);

[mT1,nT1]=size(T1);
col=T1(:,1);row=T1(1,:);
arow=col'; acol=row'; 
y=[y0;zerosq(n-1,1)];
Ty=fastQmuliply_Topelvec(acol,arow,y); 

%% ------------------TEST-------------------

tol = 1e-7;  maxit = 5000;
 E=eye(n); ecol=E(:,1);


 
%=========== without precondtiner but fast compuation ==========

e1=zerosq(nT1,1);e1(1)=quaternion(1,0,0,0);
cT1=fastQmuliply_Topelvec(col,row,e1); 

Tcol=fastQmuliply_Topelvec(acol,arow,cT1)./m;
Trow=Tcol';

zov=zeros(n,1); Ecol=quaternion(ecol,zov,zov,zov);
 
t1=tic;
[xend_fast,nomres_end_fast,iter_fast ]=fastQuaPCG(Tcol,Trow,Ty,Ecol,tol,maxit);
time_fast=toc(t1)

 
nomres_end_fast
iter_fast
 
% =========== with preconditoner and fast compuation =========

rvec=Tcol; tcn= [zerosq(d+1,1);X(1:m-1)]; rn=col'*tcn;
rvec_inv=[Trow(2:n),rn];rvec_inv=rvec_inv.';
Lrvec=[rvec;rvec_inv];

c=TCpreconder(Lrvec);

t2=tic;
[xend_prefast,nomres_end_prefast,iter_prefast ]=fastQuaPCG(Tcol,Trow,Ty,c,tol,maxit);
time_prefast=toc(t2)
 
 
nomres_end_prefast
iter_prefast


%norm(TT1*xend_prefast-Ty)


 
 
 