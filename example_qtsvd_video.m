%% The code tests the images from yuvfiles:
% the videos are downloaded from 
% http://trace.eas.asu.edu/yuv/ 

 
clc
clear
addpath('qtfm')
addpath('tensorlab4.0beta')
%% ======= load datasets ======

 load data_mobile


%% ====== generate tensor T from image sets ========


  m=length(T);
  M1=T{1};
  [n1,n2,d]=size(M1); 
  T4d=zeros(n1,n2,d,m);
 for i=1:m
 Mi=T{i};
 Mi=im2double(Mi);
 T4d(:,:,:,i)=Mi;
 end
 
[n1,n2,d,m]=size(T4d);
 
T0=zeros(n1,n2,m);
T1=reshape(T4d(:,:,1,:),[n1,n2,m]);
T2=reshape(T4d(:,:,2,:),[n1,n2,m]);
T3=reshape(T4d(:,:,3,:),[n1,n2,m]);
TQ0=quaternion(T0,T1,T2,T3);

 [TQ,meanTQ]=preproce(TQ0);
   

  

%% ====== Test methods ==========
  
 
[n1,n2,m]=size(TQ0);
mu=quaternion(1/sqrt(3),1/sqrt(3),1/sqrt(3));
alpha=quaternion(1/sqrt(2),-1/sqrt(2),0);
beta=mu*alpha; beta=v(beta);

axis=[mu,alpha,beta];

r=20
%% ---------- cross-quaternion T-SVD ------------

tic
[UT_ctsvd,ST_ctsvd,VT_ctsvd,Tap_ctsvd,err_ctsvd,rel_err_ctsvd,U_cinF,S_cinF,V_cinF]=tsvd_crossquten(TQ,r,axis);
toc

 Tap0_ctsvd=zerosq(n1,n2,m);
for t=1:m
    Tap0_ctsvd(:,:,t)=Tap_ctsvd(:,:,t)+meanTQ;
end
   [rel_err0_ctsvd,err0_ctsvd]=normQ3d(TQ0,Tap0_ctsvd);