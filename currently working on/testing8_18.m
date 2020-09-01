clear



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing hamitonian from miaojianjian PRL using Sz Sx %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%this is the main function to call other functions to do dmrg sweeps
%the parameters are set in this page  
%system parameter
format long
t = 4;
%perturbation
%tt = t-0.0001;
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
D = [0,1;1,0];%Sx
Sz = [1,0;0,-1];%Sz
%MPO formulation of Hamiltonian
%Hamiltonian tensor for each chain site
Hloc = zeros(4,4,2,2);
Hloc(1,1,:,:) = eye(2);
Hloc(2,1,:,:) = D;
Hloc(3,1,:,:) = Sz;
Hloc(4,1,:,:) = mu/2 * Sz;
Hloc(4,2,:,:) = -4*t * D;
Hloc(4,3,:,:) = 4*U * Sz;
Hloc(4,4,:,:) = eye(2);
%pertubation
Hloc_p = Hloc;
Hloc_p(4,2,:,:) = -4*tt * D;
Hloc_p = permute(Hloc_p,[1 3 2 4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hloc = permute(Hloc,[1 3 2 4]); % leg order: left-bottom-right-top%%%%(DONT CHANGE THIS)%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hs = cell(1,N);
Hs(1:N) = {Hloc};
%perturbation
Hs(2:2:N) = {Hloc_p};
% choose the last components of the left leg
Hs{1} = Hs{1}(end,:,:,:);
% choose the first components of the right leg
Hs{end} = Hs{end}(:,:,1,:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[M0,E0,Eiter0] = DMRG_1site (Hs,Nkeep,Nsweep);
%[A,sWeight,B,Eiter] = My_dmrg_2site (Hs,Nkeep,Nsweep,OPTS);
[A,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %compute correlation function for the nearest-neighbour spins
 SS = zeros(N,N);Mark = [];
D2 = 2*[0,-1;1,0]; D = D*2;
for itN = (1:N-1)%<i*lambda1_i*lambda2_j> when i>j
    for itN2 = ((itN+1:N))
        T= ncon({A{itN},D,Z,A{itN}},{[1,2,-1],[3,4],[2,4],[1,3,-2]});
        if itN2-itN > 1
            for itN3 = (itN+1:itN2-1)
                %T = ncon({T,A{itN3},A{itN3}},{[1,2],[1,3,-1],[2,3,-2]});
                T = ncon({T,A{itN3},Z,A{itN3}},{[1,2],[1,3,-1],[4,3],[2,4,-2]});
            end
        end
        if itN2 == N
            T = ncon({T,A{itN2},D2,A{itN2}},{[1,2],[1,3],[4,3],[2,4]});
        else
            T = ncon({T,A{itN2},D2,A{itN2}},{[1,2],[1,3,-1],[4,3],[2,4,-2]});
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


for itN = (2:N)
 T= ncon({A{itN-1},Sz,A{itN-1}},{[1,2,-1],[3,2],[1,3,-2]});
 for itN2 = (itN:N)
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




