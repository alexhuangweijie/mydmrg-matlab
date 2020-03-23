


clear
j.t=10;
j.delta=0;
j.mu=0;
j.u=0;

Psi = [1;5;3;14];
BLOCK1.basis_size = 2;
BLOCK1.I = eye(2);
BLOCK1.C = [0,0;0,2];%lambda1
BLOCK1.Cdag = [0,0;2,0];%lambda2
BLOCK1.H = [0,0;0,0];
BLOCK1.N = [0,0;0,2];
BLOCK2.basis_size = 2;
BLOCK2.I = eye(2);
BLOCK2.C = [0,0;1,0];%lambda1
BLOCK2.Cdag = [0,1;0,0];%lambda2
BLOCK2.H = [0,0;0,0];
BLOCK2.N = [0,0;0,1];

[H1,H4,H_super] = BlockBuilding_hubbard(BLOCK2,BLOCK2,BLOCK2,BLOCK2,j,0);
 H = H_super*H_super.';
% [V1,V2,c2,q] = lanczos(H);
%[V, V3,c1,q] = svd_R(H1,H4,j,0);
% [a1,a2] =  eig(H);
A = H;

% EYE=q'*q;
%T=q'*A*q;
%eigval1 = eigs(T,1,'SR');
%[V2,eigval2] = eigs(A,1,'SR');
%delta = eigval2-eigval1
% delta2 = eigval2-c1
% det(EYE)
% toc


% c = mykron(BLOCK1,BLOCK2,Psi,j)
% d  =  kron(BLOCK1.H,BLOCK2.I) + kron(BLOCK1.I,BLOCK2.H) ...
%         +(0.5)*(j.t+j.delta)*kron(BLOCK1.Cdag,BLOCK2.C)...
%         -(0.5)*(j.t-j.delta)*kron(BLOCK1.C,BLOCK2.Cdag)...
%         -j.u * kron(BLOCK1.N,BLOCK2.N) ...
%         +(0.5)*j.mu*(kron(BLOCK1.N,BLOCK2.I)+kron(BLOCK1.I,BLOCK2.N));
% e = d*d.'*Psi
% f = d.'*d*Psi
% 
% function [e] = mykron(a,b,p,J)
% d = reshape(p,b.basis_size,a.basis_size);
% d = d.';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %pay attention   kron(a,b)*p is beccomes a*d*b.'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h1 = a.H - 0.5 * J.mu * a.N;
% h2 = b.H - 0.5 * J.mu * b.N;
% j1 = (0.5)*(J.t+J.delta);
% j2 = -(0.5)*(J.t-J.delta);
% j3 = (0.5/1i)*J.mu;
% e =   h1.'*h1*d + h1.'*d*(h2).' + j1*h1.'*a.Cdag*d*(b.C).' + j2*h1.'*a.C*d*(b.Cdag).' + j3*h1.'*a.N*d*(b.N).'...
%     +h1*d*(h2).' + d*(h2.'*h2).'+j1*a.Cdag*d*(h2.'*b.C).'+j2*a.C*d*(h2.'*b.Cdag).'+j3*a.N*d*(h2.'*b.N).'...
%     +j1*( (a.Cdag).'*h1*d*((b.C).').' + (a.Cdag).'*d*((b.C).'*h2).'+j1*(a.Cdag).'*a.Cdag*d*((b.C).'*b.C).'...
%     +j2*(a.Cdag).'*a.C*d*((b.C).'*b.Cdag).'+j3*(a.Cdag).'*a.N*d*((b.C).'*b.N).')...
%     +j2*( (a.C).'*h1*d*((b.Cdag).').' + (a.C).'*d*((b.Cdag).'*h2).'+j1*(a.C).'*a.Cdag*d*((b.Cdag).'*b.C).'...
%     +j2*(a.C).'*a.C*d*((b.Cdag).'*b.Cdag).'+j3*(a.C).'*a.N*d*((b.Cdag).'*b.N).')...
%     +j3*((a.N).'*h1*d*((b.N).').' + (a.N).'*d*((b.N).'*h2).'+j1*(a.N).'*a.Cdag*d*((b.N).'*b.C).'...
%     +j2*(a.N).'*a.C*d*((b.N).'*b.Cdag).'+j3*(a.N).'*a.N*d*((b.N).'*b.N).');
% e = e.';
% e = e(:);
% end
%function [V1,V2,c1,q]=lanczos(A)

%clear
% tic
% Create a random symmetric matrix
 

