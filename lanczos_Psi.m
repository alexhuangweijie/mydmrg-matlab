function [V,c1]=lanczos_Psi(H1,H2,H3,H4,J,Psi_lanczos)

%clear
% tic
% Create a random symmetric matrix
 format long
%  D=900; A = rand(D);A = A + A';
% d = size(A);
% D = d(1);
S1 = H1.basis_size;
S2 = H2.basis_size;
S3 = H3.basis_size;
S4 = H4.basis_size;
D = S1*S2*S3*S4;
%Iteration with j=0
% r0=rand(D,1);
r0=Psi_lanczos;
b0=sqrt(r0'*r0);
q(:,1)=r0/b0;
a(1)=q(:,1)'*trial2(H1,H2,H3,H4,J,q(:,1));
% Iteration with j=1
r(:,1)=trial2(H1,H2,H3,H4,J,q(:,1))-a(1)*q(:,1);
b(1)=sqrt(r(:,1)'*r(:,1));
q(:,2)=r(:,1)/b(1);
a(2)=q(:,2)'*trial2(H1,H2,H3,H4,J,q(:,2));
% Iteration with j>=2
j=1;b(2)=9;c0 = rand;c1 = 0;
T1 = zeros(j,j);T1(1,1)=a(1);
while (abs(c1-c0)>10^(-16))&&(j<D)
    q(:,j+1)=r(:,j)/b(j);
    j = j+1;
    a(j)=q(:,j)'*trial2(H1,H2,H3,H4,J,q(:,j));
    r(:,j)=trial2(H1,H2,H3,H4,J,q(:,j))-a(:,j)*q(:,j)-b(j-1)*q(:,j-1);
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