clear  
clc
addpath('qtfm')
addpath('tensorlab4.0beta')




%% ------- generate a random quaternion tensor --------

n1=20;n2=30;m=10;

TQ=randq(n1,n2,m);

mu=quaternion(1/sqrt(3), 1/sqrt(3),1/sqrt(3));
alpha=quaternion(-1/sqrt(2),1/sqrt(2),0);
beta=mu*alpha; beta=v(beta);

axis=[mu,alpha,beta];

%% ---------- cross-quaternion T-SVD (QT-SVD) ------------

r=20;
tic
[UT_ctsvd,ST_ctsvd,VT_ctsvd,Tap_ctsvd,err_ctsvd,rel_err_ctsvd,UinF,SinF,VinF]=tsvd_crossquten(TQ,r,axis);
toc


dis=normQ3d(Tap_ctsvd, TQ)