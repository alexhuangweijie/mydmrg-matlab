function [V,c1]=lanczos(A)

%clear
% tic
% Create a random symmetric matrix
 format long
%  D=900; A = rand(D);A = A + A';
d = size(A);
D = d(1);
%Iteration with j=0
r0=rand(D,1);
b0=sqrt(r0'*r0);
q(:,1)=r0/b0;
a(1)=q(:,1)'*A*q(:,1);
% Iteration with j=1
r(:,1)=A*q(:,1)-a(1)*q(:,1);
b(1)=sqrt(r(:,1)'*r(:,1));
q(:,2)=r(:,1)/b(1);
a(2)=q(:,2)'*A*q(:,2);
% Iteration with j>=2
j=1;b(2)=9;c0 = rand;c1 = 0;
T1 = zeros(j,j);T1(1,1)=a(1);
while (abs(c1-c0)>10^(-16))&&(j<D)
    q(:,j+1)=r(:,j)/b(j);
    j = j+1;
    a(j)=q(:,j)'*A*q(:,j);
    r(:,j)=A*q(:,j)-a(:,j)*q(:,j)-b(j-1)*q(:,j-1);
    b(j)=sqrt(r(:,j)'*r(:,j));
    T1(j,j)=a(j);
    T1(j,j-1)=b(j-1);
    T1(j-1,j)=b(j-1);
    c0 = c1;
    [V1,c1] = eig(T1);
    [ ~ , index ] = sort( diag ( c1 ) , 'descend' );
    V1 = V1(:,index(end));
    c1 = min(min(c1));
end
V = zeros(D,1);
V = q*V1;
end
% EYE=q'*q;
%T=q'*A*q;
%eigval1 = eigs(T,1,'SR');
%[V2,eigval2] = eigs(A,1,'SR');
%delta = eigval2-eigval1
% delta2 = eigval2-c1
% det(EYE)
% toc