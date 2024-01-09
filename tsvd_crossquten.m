 %% this function generates quaternion [U,S,V] which computed by our cross-block QT-SVD  
 %% Input: quaternion tensor  T =T0+ i T1+ j T2+ k T3. 

 %% Reference: J. Pan and M. Ng, "Block Diagonalization of Quaternion Circulant Matrices with Applications "   https://arxiv.org/abs/2302.04086  
 
function [UT,ST,VT,Tap,err,rel_err,U_inF,S_inF,V_inF]=tsvd_crossquten(TQ,r,axis)

[n1,n2,m]=size(TQ);
S0=s(TQ);S1=x(TQ);S2=y(TQ);S3=z(TQ);

%----- Rewrite cirulant matrix in (mu,alpha,beta)----

mu=axis(1); alpha=axis(2); beta=axis(3);
V_mu=vec_qua(mu);
V_alpha=vec_qua(alpha);
V_beta=vec_qua(beta);

A0=S0;
A1=V_mu(1)*S1+V_mu(2)*S2+V_mu(3)*S3;
A2=V_alpha(1)*S1+V_alpha(2)*S2+V_alpha(3)*S3;
A3=V_beta(1)*S1+V_beta(2)*S2+V_beta(3)*S3;
 

% ----- generate tensors diAT and fdiAT -----

diAT=qfft_ten(quaternion(A0),mu)+qfft_ten(quaternion(A1),mu)*mu;
fdiAT=qfft_ten(quaternion(A2),mu)*alpha+qfft_ten(quaternion(A3),mu)*beta;
 
%AA=qfft_ten(TQ,mu);
%dis=normQ3d(diAT+fdiAT,AA)

% ----- Q-SVD in blocks -----

%-------- store ui,sigmai vi------------------
%--------------------------------

U_inF={}; S_inF={};  V_inF={};

%---------------------------------
%--------------------------------------------

U=zerosq(n1,r,m); Sigma=zerosq(r,r,m); V=zerosq(n2,r,m);
fTap=zerosq(n1,n2,m);

mat1_qfft=diAT(:,:,1)+fdiAT(:,:,1);

%baby=norm(AA(:,:,1)-mat1_qfft)

[u1,sigma1,v1]=svd_qua(mat1_qfft,r);
U(:,:,1)=u1; Sigma(:,:,1)=sigma1; V(:,:,1)=v1;

%-----------------------------------------------------
U_inF{1}=u1; S_inF{1}=sigma1; V_inF{1}=v1;
%-----------------------------------------------------


fTap(:,:,1)=u1*sigma1*v1';

remain=rem(m, 2);
if remain==1
    num2B=(m-1)/2;      
else
    num2B=(m-2)/2;
    
matd_qfft=diAT(:,:,2+num2B)+fdiAT(:,:,2+num2B);
[ud,sigmad,vd]=svd_qua(matd_qfft,r);
U(:,:,2+num2B)=ud; Sigma(:,:,2+num2B)=sigmad; V(:,:,2+num2B)=vd;

%---------------------------------------------------------------
U_inF{2+num2B}=ud; S_inF{2+num2B}=sigmad; V_inF{2+num2B}=vd;
%----------------------------------------------------------------


fTap(:,:,2+num2B)=ud*sigmad*vd';
end



for i=1:num2B
    mati_qfft=zerosq(2*n1,2*n2); 
mati_qfft(1:n1,1:n2)= diAT(:,:,i+1);
mati_qfft(n1+1:2*n1,n2+1:2*n2)= diAT(:,:,m+1-i);
mati_qfft(1:n1,n2+1:2*n2)= fdiAT(:,:,i+1);
mati_qfft(n1+1:2*n1,1:n2)= fdiAT(:,:,m+1-i);

[uui,ssi,vvi]=svd_qua(mati_qfft,2*r);

%---------------------------------------------------------------
U_inF{i+1}=uui; S_inF{i+1}=ssi; V_inF{i+1}=vvi;
%----------------------------------------------------------------
 

%====================================================
% % ui2=uui(1:n1,1:r); show(ui2)
% % ui1=uui(1:n1,r+1:2*r); show(ui1)
% 
% ui1=vvi(n1+1:2*n1,1:r)
% ui2=vvi(n1+1:2*n1,r+1:2*r)
% 
% suil=s(ui1)
% xui1=V_mu(1)*x(ui1)+V_mu(2)*y(ui1)+V_mu(3)*z(ui1)
% yui1=V_alpha(1)*x(ui1)+V_alpha(2)*y(ui1)+V_alpha(3)*z(ui1)
% zui1=V_beta(1)*x(ui1)+V_beta(2)*y(ui1)+V_beta(3)*z(ui1)

%=================================================
U(:,:,i+1)=uui(1:n1,1:r)+uui(1:n1,r+1:2*r);
Sigma(:,:,i+1)=ssi(1:r,1:r);
V(:,:,i+1)=vvi(1:n2,1:r)+vvi(1:n2,r+1:2*r);

U(:,:,m+1-i)=uui(n1+1:2*n1,1:r)+uui(n1+1:2*n1,r+1:2*r);
Sigma(:,:,m+1-i)=ssi(r+1:2*r,r+1:2*r);
V(:,:,m+1-i)=vvi(n2+1:2*n2,1:r)+vvi(n2+1:2*n2,r+1:2*r);

 
matapp=uui*ssi*vvi';
fTap(:,:,i+1)=matapp(1:n1,1:n2)+matapp(1:n1,n2+1:2*n2);
fTap(:,:,m+1-i)=matapp(n1+1:2*n1,1:n2)+matapp(n1+1:2*n1,n2+1:2*n2);

end

% ***** generate tensor approximation of T *****
Tap=ifften(fTap,mu);

% ***** check res *****
sTap3=tens2mat(s(Tap),3);
xTap3=tens2mat(x(Tap),3);
yTap3=tens2mat(y(Tap),3);
zTap3=tens2mat(z(Tap),3);

sT3=tens2mat(S0,3);
xT3=tens2mat(S1,3);
yT3=tens2mat(S2,3);
zT3=tens2mat(S3,3);

MaxTap=quaternion(sTap3,xTap3,yTap3,zTap3);
MaxT=quaternion(sT3,xT3,yT3,zT3);

 
difT0_compress=MaxT-MaxTap;
err=norm(difT0_compress,'fro');
rel_err=err/norm(MaxT,'fro');
 

UT= ifften(U,mu);
ST= ifften(Sigma,mu);
VT= ifften(V,mu);

end

%% =========================== functions ==================================
%% function: vector of pure quaternion
 function a=vec_qua(qa)
 a=[qa.x,qa.y,qa.z];
 end
 
 %% function: quaternion fft for 3-order quaternion tensor on its 3rd dimension.
function [fT,M3_ufd]=qfft_ten(T,mu)
addpath('qtfm')
[n1,n2,m]=size(T);  
sM3=tens2mat(s(T),3);
xM3=tens2mat(x(T),3);
yM3=tens2mat(y(T),3);
zM3=tens2mat(z(T),3);

M3_ufd=quaternion(sM3,xM3,yM3,zM3);

fM3=qfft(M3_ufd, mu, 'L');  

sfT=mat2tens(s(fM3),[n1,n2,m],3); %   use fT=fft(T,[],3) for standard fft
xfT=mat2tens(x(fM3),[n1,n2,m],3); 
yfT=mat2tens(y(fM3),[n1,n2,m],3); 
zfT=mat2tens(z(fM3),[n1,n2,m],3); 
fT=quaternion(sfT,xfT,yfT,zfT);
end

%% the following function gives tensor obtained from iqfft on a tensor.
function iTen=ifften(Ten,mu)

[n1,n2,n3]=size(Ten);
sMtens=tens2mat(s(Ten),3);
xMtens=tens2mat(x(Ten),3);
yMtens=tens2mat(y(Ten),3);
zMtens=tens2mat(z(Ten),3);

Qtens=quaternion(sMtens,xMtens,yMtens,zMtens);
iqTen=iqfft(Qtens, mu, 'L');

STen=mat2tens(s(iqTen),[n1,n2,n3],3);
XTen=mat2tens(x(iqTen),[n1,n2,n3],3);
YTen=mat2tens(y(iqTen),[n1,n2,n3],3);
ZTen=mat2tens(z(iqTen),[n1,n2,n3],3);

iTen=quaternion(STen,XTen,YTen,ZTen);

end

%% The following function gives quaternion svd for a quaternion matrix

function [u,sigma,v]=svd_qua(M,r)
[u0,s0,v0]=svd(M,0); 
u=u0(:,1:r); sigma=s0(1:r,1:r); v=v0(:,1:r);
end



