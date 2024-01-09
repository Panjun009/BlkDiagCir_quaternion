%% ------- The code tests the images from yuvfiles:---------
% the videos are downloaded from 
% http://trace.eas.asu.edu/yuv/
%% ------------ it is a revised combinated version of Video_test2 and result_video, for more details, go for them--------

%clc
clear all
addpath('qtfm')
addpath('tensorlab4.0beta')

% -----  Load datasets ------
 
  load TQ_mobile

  % load TQ_suzie  
  % load TQ_news
  %load TQ_modaugh   
  % load TQ_missA   
  %   load TQ_coast 
  % load TQ_akiyo
  % load TQ_cont

 % ------- axis -------

 [n1,n2,m]=size(TQ0);
 mu=quaternion(1/sqrt(3),1/sqrt(3),1/sqrt(3));
 alpha=quaternion(1/sqrt(2),-1/sqrt(2),0);
 beta=mu*alpha; beta=v(beta);

 axis=[mu,alpha,beta];

r=10;


% ================== Test Methods =================

% ---------- cross-quaternion T-SVD ------------

tic
[UT_ctsvd,ST_ctsvd,VT_ctsvd,Tap_ctsvd,err_ctsvd,rel_err_ctsvd,U_cinF,S_cinF,V_cinF]=tsvd_crossquten(TQ,r,axis);
toc

%  Tap0_ctsvd=zerosq(n1,n2,m);
% for t=1:m
%     Tap0_ctsvd(:,:,t)=Tap_ctsvd(:,:,t)+meanTQ;
% end
% 
%  [rel_err,err]=normQ3d_pure(TQ0,Tap0_ctsvd);



 %% ----------  quaternion Tsvd ---------
mu0=quaternion(1,0,0);
tic
[UT_qsvd,ST_qsvd,VT_qsvd,Tap_qsvd,err_qsvd,rel_err_qsvd]=tsvd_quten(TQ,r,mu0);
toc

%  Tap0_qsvd=zerosq(n1,n2,m);
% for t=1:m
%     Tap0_qsvd(:,:,t)=Tap_qsvd(:,:,t)+meanTQ;
% end
% 
% 
% [rel_err0_qsvd,err0_qsvd]=normQ3d(TQ0,Tap0_qsvd);




%% --------- tsvd on 4th-order real tensor -------
tic
[UT_4th,ST_4th,VT_4th,Tap_4th,resT_4th,resTQ_4th,errTQ_4th]=tsvd_4th(TQ,r);
toc

% Tap0_4th=zerosq(n1,n2,m);
% for t=1:m
%     Tap0_4th(:,:,t)=Tap_4th(:,:,t)+meanTQ;
% end
% 
% 
% [rel_err0_4th,err0_4th]=normQ3d(TQ0,Tap0_4th);




% save USV_mobile_80 U_cinF S_cinF V_cinF  UT_qsvd ST_qsvd VT_qsvd  UT_4th ST_4th  VT_4th
%save USV_suzie_40 U_cinF S_cinF V_cinF  UT_qsvd ST_qsvd VT_qsvd  UT_4th ST_4th  VT_4th
%save USV_cont_40 U_cinF S_cinF V_cinF  UT_qsvd ST_qsvd VT_qsvd  UT_4th ST_4th  VT_4th