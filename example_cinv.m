%% the code ici aims to compute quaternion circulant inverse by fast circulant preconditioner method
%---mu=quaternion(1,0,0);alpha=quaternion(0,1,0);beta=quaternion(0,0,1);--

clc
clear all
addpath('qtfm')

 %% ========= generate circulant quaternion matrix  and compute its inverse ========
 
 m=1000;  col=randq(m,1);

%-----------inv-------------
t0=tic;
M=QuaCirculant(col);
invM0=inv(M);
time_cinv=toc(t0);

%-----------fast-inverse-------------
 t1=tic;
 invcol=fastQinverse_circ(col);
 time_fast_cinv=toc(t1);


%% ----- check the distance of M*invM and indentity matrix

IQ=quaternion(eye(m),zeros(m),zeros(m),zeros(m));


err0_cinv_L=norm(invM0*M-IQ);
err0_cinv_R=norm(M*invM0-IQ);
dis_cinv=max(err0_cinv_L,err0_cinv_R);

invM1=QuaCirculant(invcol);
err_fast_cinv_L=norm(invM1*M-IQ);
err_fast_cinv_R=norm(M*invM1-IQ);
dis_fast_cinv=max(err_fast_cinv_L,err_fast_cinv_R) 