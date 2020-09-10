clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing hamitonian from miaojianjian PRL using Sz Sx %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is the main function to call other functions to do dmrg sweeps
%the parameters are set in this page
%system parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%using the cluster to calculate the phasediagram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format long
ppp=[];
parfor  a = 1:8
    for b = 1:4
        t = 1;
        tt = t-0.1;
        %delta = t;
        b_mu = [1, 1.5, 2, 3];
        a_Ut = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6];
        U = a_Ut(a) * t;
        mu = b_mu(b) * t;
        N = 100; % number of sites in a chain
        %DMRG parameter
        Nkeep = 50; % bond dimension
        Nsweep = 4; % number of pairs of left+right sweeps
        % OPTS.maxit = 2; % iterations of Lanczos method
        % OPTS.krydim = 4; % dimension of Krylov subspace
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
        Hs{end} = Hloc1(:,:,1,:);
        %Hs{1} = Hloc1(end,:,:,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [A,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %measurement
        T= ncon({A{1},Sx,Z,A{1}},{[1,3,-2],[2,4],[4,3],[1,2,-1]});
        for itN3 = (2:N-1)
            T = ncon({T,A{itN3},Z,A{itN3}},{[1,2],[2,3,-2],[4,3],[1,4,-1]});
        end
        T = ncon({T,A{N},Sy,A{N}},{[1,2],[2,3],[4,3],[1,4]});
        ppp(a,b) = T;
    end
end
save phasediagram
quit



