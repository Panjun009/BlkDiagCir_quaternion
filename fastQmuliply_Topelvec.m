function Th=fastQmuliply_Topelvec(col,row,h)
%% THis function computes the product of a quaternion Toeplitz matrix (col is its first column, row is its first row) and  a quaternin vector (h) in an efficient way. 
%% Note Here We only consider the situation of (i,j,k). For other situation (mu, alpha,beta), one need to transform the circulant and vector to the corresponding axis.

%mu=quaternion(1,0,0);alpha=quaternion(0,1,0);beta=quaternion(0,0,1);

m=length(col); 
row1=[row(1), flip(row(2:end))];
col1=row1.';



h1=[h;zerosq(m,1)];
circol=[col;col1];

Ch=fastQmuliply_circvec(circol,h1);

Th=Ch(1:m);
