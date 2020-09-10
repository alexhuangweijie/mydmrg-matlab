clear
%testing the measurement part 
%try to use the ncon to get the same result as the updateleft funtion



% system parameter
J = -1; % coupling strength
N = 40; % number of sites in a chain
% DMRG parameter
Nkeep = 20; % bond dimension
Nsweep = 2; % number of pairs of left+right sweeps
% Local operators
[S,I] = getLocalSpace('Spin',1/2);
% % MPO formulation of Hamiltonian
% Hamiltonian tensor for each chain site
Hloc = cell(4,4);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = squeeze(S(:,1,:));
Hloc{3,1} = squeeze(S(:,2,:));
Hloc{4,2} = J*(Hloc{2,1}');
Hloc{4,3} = J*(Hloc{3,1}');
Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));
Hloc = permute(Hloc,[3 1 4 2]); % leg order: left-bottom-right-top
% full chain
Hs = cell(1,N);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(end,:,:,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,1,:); % choose the first components of the right leg
[M_GS,E_GS,Eiter] = DMRG_2site (Hs,Nkeep,Nsweep);

% compute correlation function for the nearest-neighbour spins
SS = zeros(1,N-1);
for itN = (2:N)
    T = updateLeft([],[],M_GS{itN-1},S(:,1,:),3,M_GS{itN-1});
T = updateLeft(T,3,M_GS{itN},S(:,2,:),3,M_GS{itN});
for itN2 = ((itN+1):N)
T = updateLeft(T,2,M_GS{itN2},[],[],M_GS{itN2});
end
SS(itN-1) = T;
end
SS = SS*2; % *2 due to 1/sqrt(2) in S(:,1,:) and S(:,2,:)
% exact relation
SS_exact = (((-1).^(1:N-1))./sin((2*(1:N-1)+1)*pi/2/(N+1)) - ...
1./sin(pi/2/(N+1)))/(-2*(N+1));
figure;
plot((1:N-1),SS,'-',(1:N-1),SS_exact,'--','LineWidth',1);
legend({'DMRG','Exact'});
set(gca,'FontSize',13,'LineWidth',1);
xlabel('$\ell$','Interpreter','latex');
ylabel('$\langle \hat{S}^+_\ell \hat{S}^-_{\ell+1} \rangle$', ...
'Interpreter','latex');
xlim([1 N-1]);
grid on;


SS2 = zeros(N-1,N-1);Mark = [];
S1 = squeeze(S(:,1,:));
S2 = squeeze(S(:,2,:));

for itM = 2:N
    T2 = ncon({M_GS{itM-1}, S1, M_GS{itM-1}},{[1,2,-2],[3,2],[1,3,-1]});
    T2 = ncon({T2, M_GS{itM}, S2, M_GS{itM}},{[1,2],[2,3,-2],[4,3],[1,4,-1]}); 
for itM2 = itM+1 : N  
    T2 = ncon({T2, M_GS{itM2}, M_GS{itM2}},{[1,2],[2,3,-2],[1,3,-1]});
end
SS2(itM-1) = T2;
end
plot(SS2*2);



















