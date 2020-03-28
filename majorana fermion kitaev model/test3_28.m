clear
J.t=1;
J.delta=0;
J.mu=0;
J.u=0;
D = 2;
sysBlock.Dim = D;sysBlock.basis_size = D;
sysBlock.H = rand(D,D);
sysBlock.C = rand(D,D);
sysBlock.Cdag = rand(D,D);
sysBlock.N = rand(D,D);
sysBlock.I = eye(sysBlock.Dim);
envBlock.Dim = D;envBlock.basis_size = D;
envBlock.H = rand(D,D);
envBlock.C = rand(D,D);
envBlock.Cdag = rand(D,D);
envBlock.N = rand(D,D);
envBlock.I = eye(envBlock.Dim);
BLOCK.basis_size = 2;
BLOCK.I = eye(2);
BLOCK.C = [0,0;1,0];%lambda1
BLOCK.Cdag = [0,1;0,0];%lambda2
BLOCK.H = [0,0;0,0];
BLOCK.N = [0,0;0,1];%lambda1*lambda2
EMPTY.basis_size = 1;
a = sysBlock; H1 = sysBlock;
b = envBlock; H2 = envBlock;
J1 = 0.5 * (J.t + J.delta);
J2 = -0.5 * (J.t - J.delta);
J3 = 1i * J.u;
for i = 1:4
    [H1,H2,H]=BlockBuilding_hubbard(H1,BLOCK,BLOCK,H2,J);
end
% analytic situation
H_analytic = kron(H1.H,H2.I) + kron(H1.I,H2.H) ...
    +J1 * kron(H1.Cdag,H2.C)...
    +J2 * kron(H1.C,H2.Cdag)...
    +J3 * kron(H1.N,H2.N) ...
    -(0.5)*J.mu*(kron(H1.N,H2.I)+kron(H1.I,H2.N));
%analytic solution 
[~,eigenvalue_analyticL] = eigs(H_analytic * H_analytic', 1, 'sr');
[~,eigenvalue_analyticR] = eigs(H_analytic' * H_analytic, 1, 'sr');
%calss function
Tensor = HamTensorContraction();
Tensor.J = J;
Tensor.sysBlock = H1;
Tensor.envBlock = H2;
[Psi_GSL, eigenvalue_classL] = eigs(@Tensor.tensorContractionL, ...
    length(H1.H) *length(H2.H), 1, 'sr');
[Psi_GSR, eigenvalue_classR] = eigs(@Tensor.tensorContractionR, ...
    length(H1.H) *length(H2.H), 1, 'sr');
%lanczos method function
[~, ~, eigenvalue_lanczosL] = svd_L(H1, H2, J, 0);
[~, ~, eigenvalue_lanczosR] = svd_R(H1, H2, J, 0);



%-------------------------------------------------------------------------%
%testing the mykron function
% Psi = rand(size(H_analytic,1),1);
% C_analyticL = H_analytic*H_analytic'*Psi;
% C_analyticR = H_analytic'*H_analytic*Psi;
%calss function
% Tensor = HamTensorContraction();
% Tensor.J = J;
% Tensor.sysBlock = H1;
% Tensor.envBlock = H2;
% C_classL = Tensor.tensorContractionL(Psi);
% C_classR = Tensor.tensorContractionR(Psi);
% %mykron function
% C_mykronL = mykronL(H1, H2, Psi, J);
% C_mykronR = mykronR(H1, H2, Psi, J);
% %output the error
% sum(norm(C_analyticL - C_classL))
% sum(norm(C_analyticL - C_mykronL))
% sum(norm(C_analyticR - C_classR))
% sum(norm(C_analyticR - C_mykronR))
%-------------------------------------------------------------------------%





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
    +h1*d*(h2.').' + d*(h2*h2.').'+j1*a.Cdag*d*(b.C*h2.').'+j2*a.C*d*(b.Cdag*h2.').'+j3*a.N*d*(b.N*h2.').'...
    +j1*( h1*(a.Cdag).'*d*((b.C).').' + (a.Cdag).'*d*(h2*(b.C).').'+j1*a.Cdag*(a.Cdag).'*d*(b.C*(b.C).').'...
    +j2*a.C*(a.Cdag).'*d*(b.Cdag*(b.C).').'+j3*a.N*(a.Cdag).'*d*(b.N*(b.C).').')...
    +j2*( h1*(a.C).'*d*((b.Cdag).').' + (a.C).'*d*(h2*(b.Cdag).').'+j1*a.Cdag*(a.C).'*d*(b.C*(b.Cdag).').'...
    +j2*a.C*(a.C).'*d*(b.Cdag*(b.Cdag).').'+j3*a.N*(a.C).'*d*(b.N*(b.Cdag).').')...
    +j3*(h1*(a.N).'*d*((b.N).').' + (a.N).'*d*(h2*(b.N).').'+j1*a.Cdag*(a.N).'*d*(b.C*(b.N).').'...
    +j2*a.C*(a.N).'*d*(b.Cdag*(b.N).').'+j3*a.N*(a.N).'*d*(b.N*(b.N).').');
e = e.';
e = e(:);
end
function [e] = mykronR(a,b,p,J)
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
e =  h1.'*h1*d + h1.'*d*(h2).' + j1*h1.'*a.Cdag*d*(b.C).' + j2*h1.'*a.C*d*(b.Cdag).' + j3*h1.'*a.N*d*(b.N).'...
    +h1*d*(h2.').' + d*(h2.'*h2).'+j1*a.Cdag*d*(h2.'*b.C).'+j2*a.C*d*(h2.'*b.Cdag).'+j3*a.N*d*(h2.'*b.N).'...
    +j1*( (a.Cdag).'*h1*d*((b.C).').' + (a.Cdag).'*d*((b.C).'*h2).'+j1*(a.Cdag).'*a.Cdag*d*((b.C).'*b.C).'...
    +j2*(a.Cdag).'*a.C*d*((b.C).'*b.Cdag).'+j3*(a.Cdag).'*a.N*d*((b.C).'*b.N).')...
    +j2*( (a.C).'*h1*d*((b.Cdag).').' + (a.C).'*d*((b.Cdag).'*h2).'+j1*(a.C).'*a.Cdag*d*((b.Cdag).'*b.C).'...
    +j2*(a.C).'*a.C*d*((b.Cdag).'*b.Cdag).'+j3*(a.C).'*a.N*d*((b.Cdag).'*b.N).')...
    +j3*((a.N).'*h1*d*((b.N).').' + (a.N).'*d*((b.N).'*h2).'+j1*(a.N).'*a.Cdag*d*((b.N).'*b.C).'...
    +j2*(a.N).'*a.C*d*((b.N).'*b.Cdag).'+j3*(a.N).'*a.N*d*((b.N).'*b.N).');
e = e.';
e = e(:);
end