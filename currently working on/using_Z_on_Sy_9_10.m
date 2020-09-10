clear

% Created by W.huang (September 10, 2020) : using the Z operator on Sx will lead to 
% a new operator that is not Hermitian. 
% in this file will change the method of couting Z, which means change the 
% direction of the counting the ferimions, where the Z will add on Sy and
% fill all the empty between i and j.


%system parameter
format long
t = 1;
%perturbation
tt = t-0.00001;
delta = t;
U = 0.5 * t;
mu = 1 * t;
N = 100; % number of sites in a chain
%DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 4; % number of pairs of left+right sweeps
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 4; % dimension of Krylov subspace
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
% Hloc1 = zeros(4,4,2,2);%adding perturbation
% Hloc1(1,1,:,:) = eye(2);
% Hloc1(2,1,:,:) = Sx;
% Hloc1(3,1,:,:) = Sz;
% Hloc1(4,1,:,:) = mu * Sz;
% Hloc1(4,2,:,:) = -4* tt * Sx;
% Hloc1(4,3,:,:) = 4*U * Sz;
% Hloc1(4,4,:,:) = eye(2);
% Hloc1 = permute(Hloc1,[1 3 2 4]);
% Hs(2:2:N) = {Hloc1};
% Hs{end} = Hloc1(:,:,1,:);%rotation_bond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A0,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SS = zeros(N,N);Mark = [];
 A = A0;
 % %                    1     3
%          /--------->- A ->--  -2        
%          |            | 2             
%          ^            3              
%          |            | 2                   
%          1            Sx         
%          |            | 1             
%          ^            2               
%          |            | 2              
%          \--------->- A->-- -1        
%                    1     3
for itN = (1:N-1)
    %<A| i * lambda1_i * lambda2_j |A> when i>j
    for itN2 = ((itN+1:N))
        T= ncon({A{itN},2*Sx,A{itN}},{[1,2,-1],[2,3],[1,3,-2]});       
        if itN2-itN > 1
            for itN3 = (itN+1:itN2-1)               
                T = ncon({T,A{itN3},Z,A{itN3}},{[1,2],[1,3,-1],[3,4],[2,4,-2]});
            end
        end
        if itN2 == N
            T = ncon({T,A{itN2},Z,2*Sy,A{itN2}},{[1,2],[1,3],[3,4],[4,5],[2,5]});
        else
            T = ncon({T,A{itN2},Z,2*Sy,A{itN2}},{[1,2],[1,3,-1],[3,4],[4,5],[2,5,-2]});
        end      
        if itN2 <N
            for itN4 = itN2+1:N
                T = ncon({T,A{itN4},A{itN4}},{[1,2],[1,3,-1],[2,3,-2]});
            end
        end
        SS(itN,itN2) = T;
        Mark = [Mark;itN,itN2,T];
    end
end

for itN = (2:N)
 T= ncon({A{itN-1},2*Sz,A{itN-1}},{[1,2,-1],[2,3],[1,3,-2]});
 for itN2 = (itN:N)
 T = ncon({T,A{itN2},A{itN2}},{[1,2],[1,3,-1],[2,3,-2]});
 end
 SS(itN,itN) = T;
 Mark = [Mark;itN,itN2,T];
end 

figure(2)
SS2 = abs(SS);
[X, Y] = meshgrid(1:N, 1:N);
meshz(X, Y, SS2)
zlabel('$|G_{ij}|$','Interpreter','latex','fontsize',20);
xlabel('j');
ylabel('i');
title('原始方法测量（不做U（1）变换）')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the G(I,I+1) G(I+2,I+3) G(1,L) with different phase
%
A = A0;
phase = zeros(360,7);
for itp = 1 : 360
    %     itp = 180;%set 90 degree
    % % rotation
    U_1_rotation_bond = zeros(2,2,2,2);
    %rotation angle of the U(1) rotation
    U_1_rotation_bond(1,1,2,2) = exp(-2*pi/360*itp*1i);
    U_1_rotation_bond(2,1,1,2) = exp(2*pi/360*itp*1i);
    U_1_rotation_bond(1,2,2,1) = exp(2*pi/360*itp*1i);
    U_1_rotation_bond(2,2,1,1) = exp(-2*pi/360*itp*1i);
    U_1 = U_1_rotation_bond;
    phase(itp, 1) = rotate_measurement2(A, N, U_1, 1, N);
    phase(itp, 2) = rotate_measurement2(A, N, U_1, 1,2);
    phase(itp, 3) = rotate_measurement2(A, N, U_1, 2,3);
    phase(itp, 4) = rotate_measurement2(A, N, U_1, 23,24);
    phase(itp, 5) = rotate_measurement2(A, N, U_1, 34,35);
end
figure(4)
plot(real(phase(:,1:5)))
lgd = legend('G_{1L}', 'G_{1,2}','G_{2,3}', 'G_{23,24}', 'G_{34,35}');%, 'G_{45,46}', 'G_{56,57}');
lgd.FontSize = 12;
lgd.FontWeight = 'bold';
xlabel('phase');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %compute correlation function for the nearest-neighbour spins
SS = zeros(N,N);Mark = [];
lambda2 = [0,1;1,0]; lambda1 = 2*Sx;
for itN = (1:N-1)%<i*lambda1_i*lambda2_j> when i>j
    for itN2 = ((itN+1:N))
        T= ncon({A{itN},lambda1,conj(A{itN})},{[1,2,-2],[3,2],[1,3,-1]});
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




