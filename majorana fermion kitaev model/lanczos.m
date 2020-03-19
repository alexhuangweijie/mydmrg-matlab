function [V,c1]=lanczos(H1 ,H4 ,J ,Psi)

%clear
% tic
% Create a random symmetric matrix
 format long
%  D=900; A = rand(D);A = A + A';
D = H1.basis_size*H4.basis_size;
%D = d(1);
%Iteration with j=0
if length(Psi) == D
    r0 = Psi;
else
    r0 = rand(D,1);
end
b0=sqrt(r0'*r0);
q(:,1)=r0/b0;
%a(1)=q(:,1)'*A*q(:,1);
a(1)=q(:,1)'*mykron(H1,H4,q(:,1),J);
%a(1)=q(:,1)'*A*q(:,1);
a(1)=q(:,1)'*mykron(H1,H4,q(:,1),J);
% Iteration with j=1
r(:,1)=mykron(H1,H4,q(:,1),J)-a(1)*q(:,1);
b(1)=sqrt(r(:,1)'*r(:,1));
q(:,2)=r(:,1)/b(1);
a(2)=q(:,2)'*mykron(H1,H4,q(:,2),J);
% Iteration with j>=2
j=1;b(2)=9;c0 = rand;c1 = 0;
T1 = zeros(j,j);T1(1,1)=a(1);
while (abs(c1-c0)>10^(-16))&&(j<D)
    q(:,j+1)=r(:,j)/b(j);
    j = j+1;
    a(j)=q(:,j)'*mykron(H1,H4,q(:,j),J);
    r(:,j)=mykron(H1,H4,q(:,j),J)-a(:,j)*q(:,j)-b(j-1)*q(:,j-1);
    %r(:,j)=A*q(:,j)-a(:,j)*q(:,j)-b(j-1)*q(:,j-1);
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

function [e] = mykron(a,b,p,J) 
d = reshape(p,b.basis_size,a.basis_size);
%H = kron( H1.H , H4.I ) + kron( H1.I , H4.H ) + J *( kron( H1.Cdag , H4.C ) +  kron( H1.C , H4.Cdag ));
%e = b*d*(a.');
% H1.H  = kron(H1.H,H2.I) + kron(H1.I,H2.H) ...
%         -(0.5*1i)*(J.t+j.delta)*kron(H1.Cdag,H2.C)...
%         -(0.5*1i)*(J.t-j.delta)*kron(H1.C,H2.Cdag)...
%         -J.u * kron(H1.N,H2.N) ...
%         +(0.5*1i)*J.mu*(kron(H1,N,H2.I)+kron(H1.I,H2.N));
e = d*(a.H).'+ b.H*d + -(0.5*1i)*(J.t+J.delta) * (b.C*d*(a.Cdag).') ...
    -(0.5*1i)*(J.t-J.delta) *b.Cdag*d*(a.C).'+J.u* b.N*d*(a.N).'...
    +(0.5*1i)*J.mu*(b.N*d+d*(a.N).');
e = e(:);
end
