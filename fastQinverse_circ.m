function invcol=fastQinverse_circ(col)

%% THis function computes the inverse of a quaternion circulant matrix (col is its first column) in an efficient way. We remark that its inverse is also a quaternion circulant matrix 
%% Note Here we only consider the situation of (i,j,k). For other situation (mu, alpha,beta), one needs to transform the circulant and vector to the corresponding axis.
%% Reference: J. Pan and M. Ng, "Block Diagonalization of Quaternion Circulant Matrices with Applications "   https://arxiv.org/abs/2302.04086  

mu=quaternion(1,0,0);alpha=quaternion(0,1,0);beta=quaternion(0,0,1);
m=length(col);


zov=zeros(m,1);
sc=s(col); xc=x(col); yc=y(col);zc=z(col);
 
qsc=quaternion(sc,zov,zov,zov); 
dsc0=qfft(qsc,mu,'L') ; 

qxc=quaternion(xc,zov,zov,zov); 
dxc0=qfft(qxc,mu,'L') ;
qyc=quaternion(yc,zov,zov,zov); 
dyc0=qfft(qyc,mu,'L') ;
qzc=quaternion(zc,zov,zov,zov); 
dzc0=qfft(qzc,mu,'L') ;

Termult1=dsc0+dxc0*mu;
Termult2=dyc0*alpha+dzc0*beta;

invDv=zerosq(m,1);
invAtDv=zerosq(m,1);

Bt=zerosq(2,2); invBt=zerosq(2,2); 

if rem(m, 2) == 1
%%------------ m is odd -----------------

 numblk=(m-1)*0.5+1;
% 

for t=1:numblk
    if t==1
        B1=Termult1(1)+Termult2(1);
        inB1=B1'/(B1'*B1); invDv(1)=inB1;
    else
    Bt(1,1)=Termult1(t); Bt(2,2)=Termult1(m-t+2);
    Bt(1,2)=Termult2(t); Bt(2,1)=Termult2(m-t+2);
    invBt=inv2Dim(Bt);
    invDv(t)=invBt(1,1);invDv(m-t+2)=invBt(2,2);
    invAtDv(t)=invBt(1,2); invAtDv(m-t+2)=invBt(2,1);
    end
end


else
%%---------- m is even ------------

numblk=(m-2)*0.5+2;


for t=1:numblk
    if t==1 
        B1=Termult1(1)+Termult2(1);
        inB1=B1'/(B1'*B1); invDv(1)=inB1;
    elseif t==numblk
        Bd=Termult1(numblk)+Termult2(numblk);
        inBd=Bd'/(Bd'*Bd); invDv(t)=inBd;
    else
        %tic
    Bt(1,1)=Termult1(t); Bt(2,2)=Termult1(m-t+2);
    Bt(1,2)=Termult2(t); Bt(2,1)=Termult2(m-t+2);
    invBt=inv2Dim(Bt);
    invDv(t)=invBt(1,1);invDv(m-t+2)=invBt(2,2);
    invAtDv(t)=invBt(1,2); invAtDv(m-t+2)=invBt(2,1);
    %toc
    end
end

end

%% ---------- obtain the first of the inverse of the circulant matrix -----
 a=invDv+invAtDv;
 invcol=iqfft(a,mu,'L');