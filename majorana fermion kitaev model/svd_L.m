function [V,c1]=svd_L(H1 ,H4 ,J ,Psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%svd  using lanczos method to H^T*H 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 format long
D = H1.basis_size*H4.basis_size;
%Iteration with j=0
if length(Psi) ==D
    r0 = Psi;
else
    r0 = rand(D,1);% Create a random symmetric matrix
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
while (abs(c1-c0)>10^(-8))&&(j<D)
 %while j<D
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

function [e] = mykron(a,b,p,J) %H*H^T
d = reshape(p,b.basis_size,a.basis_size);
d = d.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pay attention   kron(a,b)*p is beccomes a*d*b.' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = a.H - 0.5 * J.mu * a.N;
h2 = b.H - 0.5 * J.mu * b.N;
j1 = (0.5)*(J.t+J.delta);
j2 = -(0.5)*(J.t-J.delta);
j3 = (0.5/1i)*J.mu;
e =   h1*h1.'*d + h1.'*d*(h2).' + j1*a.Cdag*h1.'*d*(b.C).' + j2*a.C*h1.'*d*(b.Cdag).' + j3*a.N*h1.'*d*(b.N).'...
     +h1*d*(h2).' + d*(h2*h2.').'+j1*a.Cdag*d*(b.C*h2.').'+j2*a.C*d*(b.Cdag*h2.').'+j3*a.N*d*(b.N*h2.').'...
     +j1*( h1*(a.Cdag).'*d*((b.C).').' + (a.Cdag).'*d*(h2*(b.C).').'+j1*a.Cdag*(a.Cdag).'*d*(b.C*(b.C).').'...
         +j2*a.C*(a.Cdag).'*d*(b.Cdag*(b.C).').'+j3*a.N*(a.Cdag).'*d*(b.N*(b.C).').')...
     +j2*( h1*(a.C).'*d*((b.Cdag).').' + (a.C).'*d*(h2*(b.Cdag).').'+j1*a.Cdag*(a.C).'*d*(b.C*(b.Cdag).').'...
         +j2*a.C*(a.C).'*d*(b.Cdag*(b.Cdag).').'+j3*a.N*(a.C).'*d*(b.N*(b.Cdag).').')...
     +j3*(h1*(a.N).'*d*((b.N).').' + (a.N).'*d*(h2*(b.N).').'+j1*a.Cdag*(a.N).'*d*(b.c*(b.N).').'...
         +j2*a.C*(a.N).'*d*(b.Cdag*(b.N).').'+j3*a.N*(a.N).'*d*(b.N*(b.N).').');
e = e.';
e = (:);
end
