clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing hamitonian from miaojianjian PRL using Sz Sx %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding the U(1) rotation to every bonds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the parameters are set in this page
%system parameter
format long
t = 1;
tt = t;
delta = t;
U = 0.5 * t;
mu = 1 * t;
N = 100; % number of sites in a chain
%DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 4; % number of pairs of left+right sweeps
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 2; % dimension of Krylov subspace
%Local operators
Z = [-1,0;0,1];%Z
Sx = 0.5*[0,1;1,0];%Sx
Sz = 0.5*[1,0;0,-1];%Sz
Sy = 0.5*[0,1;-1,0];% i*Sy
%MPO formulation of Hamiltonian
%Hamiltonian tensor for each chain site
Hloc = zeros(4,4,2,2);
Hloc(1,1,:,:) = eye(2);
Hloc(2,1,:,:) = Sx;
Hloc(3,1,:,:) = Sz;
Hloc(4,1,:,:) = mu * Sz;
Hloc(4,2,:,:) = -4* t * Sx;
Hloc(4,3,:,:) = 4*U * Sz;
Hloc(4,4,:,:) = eye(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hloc = permute(Hloc,[1 3 2 4]); % leg order: left-bottom-right-top%%%%(DONT CHANGE THIS)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hs = cell(1,N);
Hs(1:N) = {Hloc};
Hs{1} = Hs{1}(end,:,:,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,1,:); % choose the first components of the right leg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hloc1 = zeros(4,4,2,2);
Hloc1(1,1,:,:) = eye(2);
Hloc1(2,1,:,:) = Sx;
Hloc1(3,1,:,:) = Sz;
Hloc1(4,1,:,:) = mu * Sz;
Hloc1(4,2,:,:) = -4* tt * Sx;
Hloc1(4,3,:,:) = -4*U * Sz;
Hloc1(4,4,:,:) = eye(2);
Hloc1 = permute(Hloc1,[1 3 2 4]);
Hs(2:2:N) = {Hloc1};
Hs{end} = Hloc1(:,:,1,:);rotation_bond
%Hs{1} = Hloc1(end,:,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A0,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding U1 rotation
%finding the best phase to get maximum G1L
% A = A0;
% phase = zeros(360,1);
% for itp = 1:360
% % rotation
% U_1_rotation_bond = zeros(2,2,2,2);
% %rotation angle of the U(1) rotation
% U_1_rotation_bond(1,1,1,1) = 1;
% U_1_rotation_bond(1,2,1,2) = exp(-2*pi/360*itp*1i);
% U_1_rotation_bond(2,1,2,1) = exp(-2*pi/360*itp*1i);
% U_1_rotation_bond(2,2,2,2) = 1;
% U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
% %calculate G1L
% T4 = ncon({A{1},A{2},U_1_rotation_bond,Z,2*Sx,U_1_rotation_bond,A{1},A{2}},...
%     {[1,5,6],[6,7,-4],[10,-3,7,5],[9,10],[8,9],[2,4,-2,8],[1,2,3],[3,4,-1]});
% for itp4 = 3:N-1
% T4 = ncon({T4, A{itp4}, U_1_rotation_bond, Z,U_1_rotation_bond, A{itp4}}, ...
%     {[1, 2, 3, 4], [4, 6, -4], [8, -3, 6, 3], [7, 8],  [2, 5, -2, 7], [1, 5, -1]});
% end
% T4 = ncon({T4, A{N}, U_1_rotation_bond, Z, Sy*2, U_1_rotation_bond, A{N}}, ...
%     {[1, 2, 3, 4], [4, 6, -2], [8, 10, 3, 6], [7, 8], [9, 10], [2, 5, 7, 9], [1, 5 ,-1]});
% phase(itp, 1) = T4;
% end
% plot(abs(phase))

A = A0;
phase = zeros(360,1);
for itp = 1:360
    % % rotation
    U_1_rotation_bond = zeros(2,2,2,2);
    %rotation angle of the U(1) rotation
    U_1_rotation_bond(1,1,1,1) = 1;
    U_1_rotation_bond(1,2,1,2) = exp(-2*pi/360*itp*1i);
    U_1_rotation_bond(2,1,2,1) = exp(-2*pi/360*itp*1i);
    U_1_rotation_bond(2,2,2,2) = 1;
    U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
    
    tensors_N1 = {A{1},A{2},A{1},A{2},U_1_rotation_bond,conj(U_1_rotation_bond),Z,Sz*2,};
    connects_N1 = {[7,1,5],[5,2,-4],[7,8,6],[6,9,-1],[4,-3,2,1],[10,-2,9,8],[3,4],[10,3]};
    cont_order_N1 = [3, 4, 6, 8, 9, 5, 1, 2, 7, 10];
    T4 = ncon(tensors_N1, connects_N1 ,cont_order_N1);
    for itp4 = 3:N-1
        tensors_N2 = {T4,A{itp4},A{itp4},Z,U_1_rotation_bond,conj(U_1_rotation_bond),};
        connects_N2 = {[6,7,4,5],[6,8,-1],[5,3,-4],[1,2],[2,-3,3,4],[1,-2,8,7]};
        cont_order_N2 = [2, 6, 7, 8, 4, 1, 5, 3];
        T4 = ncon(tensors_N2, connects_N2, cont_order_N2);
    end
    tensors_N3 = {T4,U_1_rotation_bond,Z,Sy*2,A{N},A{N},conj(U_1_rotation_bond)};
    connects_N3 = {[6,9,4,5],[1,2,3,4],[7,1],[8,2],[5,3],[6,10],[7,8,10,9]};
    cont_order_N3 = [8, 7, 1, 2, 3, 10, 6, 9, 4, 5];
    T4 = ncon(tensors_N3,connects_N3,cont_order_N3);
    phase(itp, 1) = T4;
end
plot(abs(phase));

