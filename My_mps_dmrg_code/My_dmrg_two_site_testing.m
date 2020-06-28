clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing   writing my own dmrg two site function%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is the main function to call other functions to do dmrg sweeps
%the parameters are set in this page  
%system parameter
format long
t = 2;
delta = t;
%G = zeros(1,5);
%for a = 1:5
U = 0.5*t;
mu = 1*t;
J_1 = 0.5*(t+delta);
J_2 = -0.5*(t-delta);
J_3 = U;
J_4 = -0.5*mu; 
N = 100; % number of sites in a chain
%DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 4; % number of pairs of left+right sweeps
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 4; % dimension of Krylov subspace
%Local operators
lambda_1 = [0,1;0,0];
lambda_2 = [0,0;1i,0];
lambda_3 = [1i,0;0,0];
%MPO formulation of Hamiltonian
%Hamiltonian tensor for each chain site
Hloc = zeros(5,5,2,2);
Hloc(1,1,:,:) = eye(2);
Hloc(2,1,:,:) = lambda_1;
Hloc(3,1,:,:) = lambda_2;
Hloc(4,1,:,:) = lambda_3;
Hloc(5,1,:,:) = lambda_3*J_4;
Hloc(5,2,:,:) = lambda_2*J_1;
Hloc(5,3,:,:) = lambda_1*J_2;
Hloc(5,4,:,:) = lambda_3*J_3;
Hloc(5,5,:,:) = eye(2);
Hloc = permute(Hloc,[1 3 2 4]); % leg order: left-bottom-right-top
%pertubation
Hs = cell(1,N);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(end,:,:,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,1,:); % choose the first components of the right leg

%[M0,E0,Eiter0] = DMRG_1site (Hs,Nkeep,Nsweep);
[A,sWeight,B,Eiter] = My_dmrg_2site (Hs,Nkeep,Nsweep,OPTS);
% %compute correlation function for the nearest-neighbour spins
SS = zeros(N,N);Mark = [];
lambda_1 = [0,0;1,0];
lambda_2 = [0,-1;0,0];
lambda_3 = [0,0;0,-1];
for itN = (1:N-1)
    for itN2 = ((itN+1:N))
        T= ncon({A{itN},lambda_1,A{itN}},{[1,2,-1],[3,2],[1,3,-2]});
        if itN2-itN > 1
            for itN3 = (itN+1:itN2-1)
                T = ncon({T,A{itN3},A{itN3}},{[1,2],[1,3,-1],[2,3,-2]});
            end
        end
        if itN2 == N
            T = ncon({T,A{itN2},lambda_2,A{itN2}},{[1,2],[1,3],[4,3],[2,4]});
        else
            T= ncon({T,A{itN2},lambda_2,A{itN2}},{[1,2],[1,3,-1],[4,3],[2,4,-2]});
        end      
        if itN2 <N
            for itN4 = itN2+1:N
                T = ncon({T,A{itN4},A{itN4}},{[1,2],[1,3,-1],[2,3,-2]});
            end
            %else
            %T = ncon({T},[1,1]);
        end
        SS(itN,itN2) = T;
        Mark = [Mark;itN,itN2,T];
    end
end
for itN = (1:N-1)
    for itN2 = ((itN+1:N))
        T= ncon({A{itN},lambda_2,A{itN}},{[1,2,-1],[3,2],[1,3,-2]});
        if itN2-itN > 1
            for itN3 = (itN+1:itN2-1)
                T = ncon({T,A{itN3},A{itN3}},{[1,2],[1,3,-1],[2,3,-2]});
            end
        end
        if itN2 == N
            T = ncon({T,A{itN2},lambda_1,A{itN2}},{[1,2],[1,3],[4,3],[2,4]});
        else
            T= ncon({T,A{itN2},lambda_1,A{itN2}},{[1,2],[1,3,-1],[4,3],[2,4,-2]});
        end      
        if itN2 <N
            for itN4 = itN2+1:N
                T = ncon({T,A{itN4},A{itN4}},{[1,2],[1,3,-1],[2,3,-2]});
            end
            %else
            %T = ncon({T},[1,1]);
        end
        SS(itN2,itN) = T;
        Mark = [Mark;itN2,itN,T];
    end
end
for itN = (2:N)
 T= ncon({A{itN-1},lambda_3,A{itN-1}},{[1,2,-1],[3,2],[1,3,-2]});
 for itN2 = ((itN):N)
 T = ncon({T,A{itN2},A{itN2}},{[1,2],[1,3,-1],[2,3,-2]});
 end
 SS(itN,itN) = T;
 Mark = [Mark;itN,itN2,T];
end 
figure(1)
SS2 = abs(SS);
[X, Y] = meshgrid(1:N, 1:N);
meshz(X, Y, SS2)
zlabel('$|G_{ij}|$','Interpreter','latex','fontsize',20);
xlabel('j');
ylabel('i');


%%
% clear   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Heisenberg model using My_dmrg_2site   
% %2020.5.26
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system parameter
% format long
% J = -1; % coupling strength
% N = 40; % number of sites in a chain
% DMRG parameter
% Nkeep = 15; % bond dimension
% Nsweep = 4; % number of pairs of left+right sweeps
% OPTS.maxit = 2; % iterations of Lanczos method
% OPTS.krydim = 4; % dimension of Krylov subspace
% Local operators
% [S,I] = getLocalSpace('Spin',1/2);
% % MPO formulation of Hamiltonian
% Hamiltonian tensor for each chain site
% Hloc = cell(4,4);
% Hloc(:) = {zeros(size(I))};
% Hloc{1,1} = I;
% Hloc{2,1} = squeeze(S(:,1,:));
% Hloc{3,1} = squeeze(S(:,2,:));
% Hloc{4,2} = J*(Hloc{2,1}');
% Hloc{4,3} = J*(Hloc{3,1}');
% Hloc{end,end} = I;
% Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));
% Hloc = permute(Hloc,[3 1 4 2]); % leg order: left-bottom-right-top
% full chain
% Hs = cell(1,N);
% Hs(:) = {Hloc};
% Hs{1} = Hs{1}(end,:,:,:); % choose the last components of the left leg
% Hs{end} = Hs{end}(:,:,1,:); % choose the first components of the right leg
%  [A,sWeight,B,Eiter] = My_dmrg_2site (Hs,Nkeep,Nsweep,OPTS);
%  measurement
%  Sp = [0,1;0,0];Sm = [0,0;1,0];
%  SS = zeros(1,N-1);
% for itN = 1:N-1
%     Atemp = ncon({A{itN},sWeight{itN+1}},{[-1,-2,1],[1,-3]});
%  T = ncon({Atemp,Sp,conj(Atemp),B{itN+1},Sm,conj(B{itN+1})},...
%      {[1,2,7],[2,3],[1,3,8],[7,5,4],[5,6],[8,6,4]});
%  SS(itN) = T;
% end
%  plot(SS)

%%
% testing the My_dmrg_2site function with Majumdar-Ghosh model
% result can match with the DMRG_goundstate_sol.pdf
% %2020.5.27
% clear
% % system parameter
% J1 = 1; % nearest-neighbour coupling strength
% J2 = 1/2; % next-nearest-neighbour coupling strength
% N = 40; % number of sites in a chain
% % DMRG parameter
% Nkeep = 2; % bond dimension
% Nsweep = 6; % number of pairs of left+right sweeps
% OPTS.maxit = 2; % iterations of Lanczos method
% OPTS.krydim = 4; % dimension of Krylov subspace
% % Local operators
% [S,I] = getLocalSpace('Spin',1/2);
% % Hamiltonian tensor for each chain site
% Hloc = cell(8,8);
% Hloc(:) = {zeros(size(I))};
% Hloc{1,1} = I;
% Hloc{2,1} = squeeze(S(:,1,:));
% Hloc{3,1} = squeeze(S(:,2,:));
% Hloc{4,1} = squeeze(S(:,3,:));
% Hloc{5,2} = I;
% Hloc{6,3} = I;
% Hloc{7,4} = I;
% Hloc{8,2} = J1*(Hloc{2,1}');
% Hloc{8,3} = J1*(Hloc{3,1}');
% Hloc{8,4} = J1*(Hloc{4,1}');
% Hloc{8,5} = J2*(Hloc{2,1}');
% Hloc{8,6} = J2*(Hloc{3,1}');
% Hloc{8,7} = J2*(Hloc{4,1}');
% Hloc{end,end} = I;
% Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));
% Hloc = permute(Hloc,[3 1 4 2]); % leg order: left-bottom-right-top
% % full chain
% Hs = cell(1,N);
% Hs(:) = {Hloc};
% Hs{1} = Hs{1}(end,:,:,:); % choose the last components of the left leg
% Hs{end} = Hs{end}(:,:,1,:); % choose the first components of the right leg
% [A,sWeight,B,Eiter] = My_dmrg_2site (Hs,Nkeep,Nsweep,OPTS);
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% ketaev model testing %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Set simulation options
% chi = 16; % maximum bond dimension
% Nsites = 100; % number of lattice sites
% OPTS.numsweeps = 4; % number of DMRG sweeps
% OPTS.display = 2; % level of output display
% OPTS.updateon = 1; % update MPS tensors
% OPTS.maxit = 20; % iterations of Lanczos method
% OPTS.krydim = 20; % dimension of Krylov subspace
% 
% %%%% Define Hamiltonian MPO (majorana fermions(kitaev model))
% chid = 2;
% t = 1;
% delta = 1;
% U = 0.5;
% mu = 1;
% J_2 = 0.5i*(t+delta);
% J_1 = -0.5i*(t-delta);
% J_3 = -U;
% J_4 = -0.5i*mu; 
% % Local operators
% lambda_1 = [0,0;1,0];
% lambda_2 = [0,1i;0,0];
% lambda_3 = [0,0;0,1i];
%  % MPO formulation of Hamiltonian
% % Hamiltonian tensor for each chain site
% Hloc = zeros(5,5,2,2);
% Hloc(1,1,:,:) = eye(2);
% Hloc(1,2,:,:) = lambda_1;
% Hloc(1,3,:,:) = lambda_2;
% Hloc(1,4,:,:) = lambda_3;
% Hloc(1,5,:,:) = lambda_3*J_4;
% Hloc(2,5,:,:) = lambda_2*J_1;
% Hloc(3,5,:,:) = lambda_1*J_2;
% Hloc(4,5,:,:) = lambda_3*J_3;
% Hloc(5,5,:,:) = eye(2);
% % %Hloc = permute(Hloc,[2 4 1 3]); % leg order: left-bottom-right-top
%  M = Hloc;
%  ML = reshape([1;0;0;0;0],[5,1,1]); %left MPO boundary
%  MR = reshape([0;0;0;0;1],[5,1,1]); %right MPO boundary
% 
% %%%% Initialize MPS tensors
% A = {};
% A{1} = rand(1,chid,min(chi,chid));
% for k = 2:Nsites
%     A{k} = rand(size(A{k-1},3),chid,min(min(chi,size(A{k-1},3)*chid),chid^(Nsites-k)));
% end
% 
% %%%% Do DMRG sweeps 
% [A,~,~,Ekeep1] = doDMRG_MPO(A,ML,M,MR,chi,OPTS);
% 
% %%%% Increase bond dim and reconverge 
% chi = 32;
% [A,sWeight,B,Ekeep2] = doDMRG_MPO(A,ML,M,MR,chi,OPTS);
% 
% %%%% Compute the correlation function , local density
% rhotwo = {};
% hamloc = reshape(kron(sP,sM) + kron(sM,sP),[2,2,2,2]);
% for k = 1:Nsites-1
%     rhotwo{k} = ncon({A{k},conj(A{k}),A{k+1},conj(A{k+1}),sWeight{k+2},sWeight{k+2}},...
%         {[1,-3,2],[1,-1,3],[2,-4,4],[3,-2,5],[4,6],[5,6]});
%     Enloc(k) = ncon({hamloc,rhotwo{k}},{[1:4],[1:4]});
% end
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% majumdar-ghosh model testing %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %结果并不正确
% %%%%% Set simulation options
% chi = 16; % maximum bond dimension
% Nsites = 100; % number of lattice sites
% OPTS.numsweeps = 4; % number of DMRG sweeps
% OPTS.display = 2; % level of output display
% OPTS.updateon = 1; % update MPS tensors
% OPTS.maxit = 30; % iterations of Lanczos method  
% OPTS.krydim = 4; % dimension of Krylov subspace 保证多少个向量是正交的
% E0 = -3/8*Nsites
% %%%% Define Hamiltonian MPO (majorana fermions(kitaev model))
% chid = 2;
% % Local operators
% lambda_z = [0.5,0;0,-0.5];
% lambda_p = [0,1;0,0];
% lambda_m = [0,0;0,1];
%  % MPO formulation of Hamiltonian
% % Hamiltonian tensor for each chain site
% Hloc = zeros(8,8,2,2);
% Hloc(1,1,:,:) = eye(2);
% Hloc(1,2,:,:) = lambda_z;
% Hloc(1,3,:,:) = lambda_m;
% Hloc(1,4,:,:) = lambda_p;
% Hloc(2,5,:,:) = eye(2);
% Hloc(3,6,:,:) = eye(2);
% Hloc(4,7,:,:) = eye(2);
% Hloc(2,8,:,:) = lambda_z;
% Hloc(3,8,:,:) = lambda_p;
% Hloc(4,8,:,:) = lambda_m;
% Hloc(5,8,:,:) = lambda_z*0.5;
% Hloc(6,8,:,:) = lambda_p*0.5;
% Hloc(7,8,:,:) = lambda_m*0.5;
% Hloc(8,8,:,:) = eye(2);
% % %Hloc = permute(Hloc,[2 4 1 3]); % leg order: left-bottom-right-top
%  M = Hloc;
%  ML = reshape([1;0;0;0;0;0;0;0],[8,1,1]); %left MPO boundary
%  MR = reshape([0;0;0;0;0;0;0;1],[8,1,1]); %right MPO boundary
% 
% %%%% Initialize MPS tensors
% A = {};
% A{1} = rand(1,chid,min(chi,chid));
% for k = 2:Nsites
%     A{k} = rand(size(A{k-1},3),chid,min(min(chi,size(A{k-1},3)*chid),chid^(Nsites-k)));
% end
% 
% %%%% Do DMRG sweeps 
% [A,~,~,Ekeep1] = doDMRG_MPO(A,ML,M,MR,chi,OPTS);
% 
% %%%% Increase bond dim and reconverge 
% chi = 32;
% [A,sWeight,B,Ekeep2] = doDMRG_MPO(A,ML,M,MR,chi,OPTS);
% 
% %%%% Compute the correlation function , local density
% rhotwo = {};
% hamloc = reshape(kron(lambda_m,lambda_p) + kron(lambda_m,lambda_p),[2,2,2,2]);
% for k = 1:Nsites-1
%     rhotwo{k} = ncon({A{k},conj(A{k}),A{k+1},conj(A{k+1}),sWeight{k+2},sWeight{k+2}},...
%         {[1,-3,2],[1,-1,3],[2,-4,4],[3,-2,5],[4,6],[5,6]});
%     Enloc(k) = ncon({hamloc,rhotwo{k}},{[1:4],[1:4]});
% end
% plot(Enloc)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %testing the right orthognality of the MPS
% T = ncon({A{50},A{50}},{[-1,2,1],[-2,2,1]})
% for i = 49:-1:20
%     T = ncon({A{i},T,A{i}},{[-1,3,1],[1,2],[-2,3,2]});
% end
% T
% %testing the right orthognality of the MPS
% T1 = ncon({A{1},A{1}},{[1,2,-1],[1,2,-2]});
% for i = 2:20
%     T1 = ncon({A{i},T1,A{i}},{[1,3,-1],[1,2],[2,3,-2]});
% end
% T1
%%