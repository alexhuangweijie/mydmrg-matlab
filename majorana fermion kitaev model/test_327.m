% By Huangweijie @ SCUT, 2020.03.027
% Running on MATLAB R2019b
% Testing the eigs funtion 
clear
J.t=1;
J.delta=0;
J.mu=0;
J.u=0;
D = 8;
% Construct Spin operator
sysBlock.Dim = D;sysBlock.basis_size = D;
sysBlock.H = eye(D,D);
sysBlock.C = eye(D,D);
sysBlock.Cdag = eye(D,D);
sysBlock.N = eye(D,D);
sysBlock.I = eye(sysBlock.Dim);

envBlock.Dim = D;
envBlock.H = eye(D,D);
envBlock.C = eye(D,D);
envBlock.Cdag = eye(D,D);
envBlock.N = eye(D,D); 
envBlock.I = eye(envBlock.Dim);
envBlock.basis_size = D;
    % Implement the tensor contraction instead of constructing superblock
    Tensor = HamTensorContraction();
    Tensor.J = J;
    Tensor.sysBlock = sysBlock;
    Tensor.envBlock = envBlock;
    
H1 = kron(sysBlock.H,envBlock.I) + kron(sysBlock.I,envBlock.H) ...
        +(0.5)*(J.t+J.delta)*kron(sysBlock.Cdag,envBlock.C)...
        -(0.5)*(J.t-J.delta)*kron(sysBlock.C,envBlock.Cdag)...
        -J.u *(1i)* kron(sysBlock.N,envBlock.N) ...
        -(0.5)*J.mu*(kron(sysBlock.N,envBlock.I)+kron(sysBlock.I,envBlock.N));
Hcal = H1'*H1; 
PSI = rand(size(Hcal,1),1);
A_analytic = Hcal * PSI;
A_tensor = Tensor.tensorContractionL(PSI);
A_mykron = mykronL(sysBlock,envBlock,PSI,J);
sum(norm(A_analytic - A_mykron))
sum(norm(A_mykron - A_tensor))
sum(norm(A_analytic - A_tensor))
%--------------------------------------------------------------------------
% [Psi_func, E_func] = eigs(@Tensor.tensorContractionL, ...
%         length(sysBlock.H) *length(envBlock.H), 1, 'sr');
% [Psi_analytic,E_analytic] = eigs(Hcal,1,'sr');
% [Psi_lanczos,E_lanczos] = svd_L(sysBlock, envBlock, J, 0);%

function [e] = mykronL(a,b,p,J) %H*H^T
d = reshape(p,b.basis_size,a.basis_size);
d = d.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pay attention   kron(a,b)*p is beccomes a*d*b.' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = a.H - 0.5 * J.mu * a.N;
h2 = b.H - 0.5 * J.mu * b.N;
j1 = (0.5)*(J.t+J.delta);
j2 = -(0.5)*(J.t-J.delta);
j3 = (1i)*J.u;
e =   h1*h1.'*d + h1.'*d*(h2).' + j1*a.Cdag*h1.'*d*(b.C).' + j2*a.C*h1.'*d*(b.Cdag).' + j3*a.N*h1.'*d*(b.N).'...
     +h1*d*(h2).' + d*(h2*h2.').'+j1*a.Cdag*d*(b.C*h2.').'+j2*a.C*d*(b.Cdag*h2.').'+j3*a.N*d*(b.N*h2.').'...
     +j1*( h1*(a.Cdag).'*d*((b.C).').' + (a.Cdag).'*d*(h2*(b.C).').'+j1*a.Cdag*(a.Cdag).'*d*(b.C*(b.C).').'...
         +j2*a.C*(a.Cdag).'*d*(b.Cdag*(b.C).').'+j3*a.N*(a.Cdag).'*d*(b.N*(b.C).').')...
     +j2*( h1*(a.C).'*d*((b.Cdag).').' + (a.C).'*d*(h2*(b.Cdag).').'+j1*a.Cdag*(a.C).'*d*(b.C*(b.Cdag).').'...
         +j2*a.C*(a.C).'*d*(b.Cdag*(b.Cdag).').'+j3*a.N*(a.C).'*d*(b.N*(b.Cdag).').')...
     +j3*(h1*(a.N).'*d*((b.N).').' + (a.N).'*d*(h2*(b.N).').'+j1*a.Cdag*(a.N).'*d*(b.C*(b.N).').'...
         +j2*a.C*(a.N).'*d*(b.Cdag*(b.N).').'+j3*a.N*(a.N).'*d*(b.N*(b.N).').');
e = e.';
e = e(:);
end


    