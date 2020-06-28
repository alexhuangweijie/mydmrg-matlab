clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing kitaev chain with jordan-wigner transform%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is the main function to call other functions to do dmrg sweeps
%the parameters are set in this page  
%system parameter
format long
t = 2;
delta = t;
U = 0.5*t;
mu = 1*t;
N = 40; % number of sites in a chain
%DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 4; % number of pairs of left+right sweeps
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 4; % dimension of Krylov subspace
%Local operators
C = [0,0;1,0];
Cdag = [0,1;0,0];
n = [1,0;0,0];
n1 = 2*n-eye(2);
n2 = n-1/2*eye(2);
%MPO formulation of Hamiltonian
%Hamiltonian tensor for each chain site
Hloc = zeros(5,5,2,2);
Hloc(1,1,:,:) = eye(2);
Hloc(2,1,:,:) = C;
Hloc(3,1,:,:) = Cdag;
Hloc(4,1,:,:) = n1;
Hloc(5,1,:,:) = mu*n2;
Hloc(5,2,:,:) = delta*Cdag - t*C;
Hloc(5,3,:,:) = delta*C + t*Cdag;
Hloc(5,4,:,:) = U*n1;
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
lambda_1 = [0,1;1,0];
lambda_2 = [0,1;-1,0];
lambda_3 = [-1,0;0,1];
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