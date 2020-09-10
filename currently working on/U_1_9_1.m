clear


% calculate the G_{i,j} by using the U_1 rotation. 
% this is the lattest edition. hwj/2020/9/01.


%system parameter
format long
t = 1;
%perturbation
tt = t-0.00001;
delta = t;
U = 1 * t;
mu = 0.5 * t;
N = 100; % number of sites in a chain
%DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 2; % number of pairs of left+right sweeps
% OPTS.maxit = 2; % iterations of Lanczos method
% OPTS.krydim = 4; % dimension of Krylov subspace
%Local operators
Z = [-1,0;0,1];%Z
Sz = 0.5*[1,0;0,-1];%Sz
Sx = 0.5*[0,1;1,0];%Sx
Sy = 0.5*[0,1;-1,0];% i*Sy
%MPO formulation of Hamiltonian
%Hamiltonian tensor for each chain site
Hloc = zeros(4,4,2,2);
Hloc(1,1,:,:) = eye(2);
Hloc(2,1,:,:) = Sx;
Hloc(3,1,:,:) = Sz;
Hloc(4,1,:,:) = mu * Sz;
Hloc(4,2,:,:) = -4*t * Sx;
Hloc(4,3,:,:) = 4*U * Sz;
Hloc(4,4,:,:) = eye(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hloc = permute(Hloc,[1 3 2 4]); % leg order: left-bottom-right-top%%%%(DONT CHANGE THIS)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hs = cell(1,N);
Hs(1:N) = {Hloc};
Hs{1} = Hs{1}(end,:,:,:);% choose the last components of the left leg
Hs{end} = Hs{end}(:,:,1,:);% choose the first components of the right leg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hloc1 = zeros(4,4,2,2);%adding perturbation
Hloc1(1,1,:,:) = eye(2);
Hloc1(2,1,:,:) = Sx;
Hloc1(3,1,:,:) = Sz;
Hloc1(4,1,:,:) = mu * Sz;
Hloc1(4,2,:,:) = -4* tt * Sx;
Hloc1(4,3,:,:) = 4*U * Sz;
Hloc1(4,4,:,:) = eye(2);
Hloc1 = permute(Hloc1,[1 3 2 4]);
Hs(2:2:N) = {Hloc1};
Hs{end} = Hloc1(:,:,1,:);%rotation_bond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A0,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the G(I,I+1) G(I+2,I+3) G(1,L) with different phase
% 
A = A0;
phase = zeros(360,3); 
for itp = 1 : 360
%     itp = 180;%set 90 degree
    % % rotation
    U_1_rotation_bond = zeros(2,2,2,2);
    %rotation angle of the U(1) rotation
    U_1_rotation_bond(1,1,1,1) = 1;
    U_1_rotation_bond(1,2,1,2) = exp(2*pi/360*itp*1i);
    U_1_rotation_bond(2,1,2,1) = exp(2*pi/360*itp*1i);
    U_1_rotation_bond(2,2,2,2) = 1;
    U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
    U_1 = U_1_rotation_bond;
    phase(itp, 1) = rotate_measurement(A, N, U_1, 1, N);
    phase(itp, 2) = rotate_measurement(A, N, U_1, 30, 31);
    phase(itp, 3) = rotate_measurement(A, N, U_1, 31, 32);
end
figure(1)
plot(real(phase))
lgd = legend('G_{1L}', 'G_{30,31}','G_{31,32}');
lgd.FontSize = 12;
lgd.FontWeight = 'bold';
xlabel('phase');


G = zeros(N,N);
%     itp = 180;%set 90 degree
% % rotation
U_1_rotation_bond = zeros(2,2,2,2);
%rotation angle of the U(1) rotation
U_1_rotation_bond(1,1,1,1) = 1;
U_1_rotation_bond(1,2,1,2) = exp(2*pi/360*90*1i);
U_1_rotation_bond(2,1,2,1) = exp(2*pi/360*90*1i);
U_1_rotation_bond(2,2,2,2) = 1;
U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
U_1 = U_1_rotation_bond;
for  index_1 = 1:N-1
    for index_2 = index_1+1:N
        G(index_1,index_2) = rotate_measurement(A, N, U_1, index_1, index_2);
    end
end

figure(2)
[X, Y] = meshgrid(1:N, 1:N);
meshz(X, Y, abs(G))
zlabel('$|G_{ij}|$','Interpreter','latex','fontsize',20);
xlabel('j');
ylabel('i');

