clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing hamitonian from miaojianjian PRL using Sz Sx %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is the main function to call other functions to do dmrg sweeps




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%adding the perturbation to he even sites
%testing the first excited state


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%the parameters are set in this page  
%system parameter
format long
t = 1;
tt = t;
delta = t;
U = 0.5 * t;
mu = 0 * t;
N = 60; % number of sites in a chain
%DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 4; % number of pairs of left+right sweeps
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 2; % dimension of Krylov subspace
%Local operators
Z = [-1,0;0,1];%Z
Sx = 0.5*[0,1;1,0];%Sx
Sz = 0.5*[1,0;0,-1];%Sz
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
[A0,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);

%first excited state
% [A1,E1,Eiter1] = DMRG_1site (Hs,Nkeep,Nsweep,A);
% A_temp = A1; 
% A1 = A;
% A = A_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding U1 rotation
A = A0;
U_1_rotation = exp(rand*1i);
U_1 = [U_1_rotation, 0; 0, 1;];
for itM = 1:N
   A{itM} = ncon({A{itM},U_1},{[-1, 1, -3],[-2, 1]});
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %compute correlation function for the nearest-neighbour spins
SS = zeros(N,N);Mark = [];
lambda2 = [0,1;-1,0]; lambda1 = 2*Sx;
for itN = (1:N-1)%<i*lambda1_i*lambda2_j> when i>j
    for itN2 = ((itN+1:N))
        T= ncon({A{itN},lambda1,Z,conj(A{itN})},{[1,2,-2],[3,4],[4,2],[1,3,-1]});
        if itN2-itN > 1
            for itN3 = (itN+1:itN2-1)
                T = ncon({T,A{itN3},Z,conj(A{itN3})},{[1,2],[2,4,-2],[3,4],[1,3,-1]});
            end
        end
        if itN2 == N
            T = ncon({T,A{itN2},lambda2,conj(A{itN2})},{[1,2],[2,4],[3,4],[1,3]});
        else
            T = ncon({T,A{itN2},lambda2,conj(A{itN2})},{[1,2],[2,4,-2],[3,4],[1,3,-1]});
        end      
        if itN2 <N
            for itN4 = itN2+1:N
                T = ncon({T,A{itN4},conj(A{itN4})},{[1,2],[2,3,-2],[1,3,-1]});
            end
        end
        SS(itN,itN2) = T;
        Mark = [Mark;itN,itN2,T];
    end
end

for itN = (2:N)
 T= ncon({A{itN-1},Sz*2,conj(A{itN-1})},{[1,3,-2],[2,3],[1,2,-1]});
 for itN2 = (itN:N)
 T = ncon({T,A{itN2},conj(A{itN2})},{[1,2],[2,3,-2],[1,3,-1]});
 end
 SS(itN-1,itN-1) = T;
 Mark = [Mark;itN,itN,T];
end 
%SS(N,N) = ncon({A{N},Sz*2,A{N}},{[1,3],[2,3],[1,2]});
figure(1)
SS2 = abs(SS);
[X, Y] = meshgrid(1:N, 1:N);
meshz(X, Y, SS2)
zlabel('$|G_{ij}|$','Interpreter','latex','fontsize',20);
xlabel('j');
ylabel('i');
figure(2)
plot(diag(SS2),'-o')

% save correlationfunction
% quit



