 function invQ=inv2Dim(Q)
%Q=2*randq(2);
% This function gives the closed form of inverse of 2*2 quaternion matrix.

invQ=zerosq(2,2);
q11=Q(1,1);q12=Q(1,2);q21=Q(2,1); q22=Q(2,2);

n11=q11*q11'; n12=q12*q12'; n21=q21*q21'; n22=q22*q22';

trq=q12'*q11*q21'*q22; nq= 2*s(trq);

detm=n11*n22+n12*n21-nq;
invQ(1,1)=n22*q11'-q21'*q22*q12';
invQ(1,2)=n12*q21'-q11'*q12*q22'; 
invQ(2,1)=n21*q12'-q22'*q21*q11'; 
invQ(2,2)=n11*q22'-q12'*q11*q21';
invQ=invQ/detm;


