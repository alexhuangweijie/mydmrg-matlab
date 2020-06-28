function [A,sWeight,B,Eiter] = My_dmrg_2site(Hs, Nkeep, Nsweep,OPTS)
N = numel(Hs);

fprintf('ground state\n');
fprintf(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2\n']);
% % % initialize with random MPS
M = cell(1,N);
for itN = (1:N)
     if itN == 1%physical_bond_size = size(Hs{itN},4);
         M{itN} = rand(1,size(Hs{itN},4),min(size(Hs{itN},4),Nkeep)); % right leg is dummy
     elseif itN == N
         M{itN} = rand(min(size(Hs{itN},4),Nkeep),size(Hs{itN},4),1); % left leg is dummy
     else
        M{itN} = rand(size(M{itN-1},3),size(Hs{itN},4),...
            min(min(Nkeep,size(M{itN-1},3)*size(Hs{itN},4)),size(Hs{itN},4)^(N-itN)));
    end
end

A = M;B = M;

% ground-state energy for each iteration
Eiter = zeros(N,2*Nsweep);
Hl = cell(1,N+2);
Hr = cell(1,N+2); 

for itN = 1:N-1                            % bring into left-canonical form
    chil = size(A{itN},1); chir = size(A{itN},3);chid = size(A{itN},2);
    [qtemp,rtemp] = qr(reshape(A{itN},[chil*chid,chir]),0);
    A{itN} = reshape(qtemp,[chil,chid,chir]);
    A{itN+1} = ncon({rtemp,A{itN+1}},{[-1,1],[1,-2,-3]})/norm(rtemp(:));
    if itN == 1
        % "remove" the left leg (1st leg) which is dummy by permuting to the last
        H2 = permute(Hs{itN},[2 3 4 1]);
        Hl{itN+1} = ncon({A{itN},H2,conj(A{itN})},{[1,3,-3],[2,-2,3],[1,2,-1]},[1,2,3]);
    else
        Hl{itN+1} = ncon({Hl{itN},A{itN},Hs{itN},conj(A{itN})},{[3,4,5],[5,2,-3],[4,1,-2,2],[3,1,-1]});
    end
end
chil = size(A{N},1); chir = size(A{N},3);chid = size(A{N},2);
[qtemp,stemp] = qr(reshape(A{N},[chil*chid,chir]),0);
A{N} = reshape(qtemp,[chil,chid,chir]);
sWeight{N +1} = stemp./sqrt(trace(stemp*stemp'));

Hl{1} = 1;Hl{end} = 1;
Hr{1} = 1;Hr{end} = 1;
disptime(['Initialize with random MPS. Energy = ',sprintf('%.7g',Hl{N+1})]);

%dmrg sweep
for itS = (1:Nsweep)
    % right -> left
    for itN = N-1:-1:1
        chil = size(A{itN},1); chir = size(A{itN+1},3); chid = size(A{itN},2);
        psiGround = reshape(ncon({A{itN},A{itN+1},sWeight{itN+2}},{[-1,-2,1],[1,-3,2],[2,-4]}),[chil*chid^2*chir,1]);
        %[Aeff, Eeff] = My_eigs_2site(psiGround,OPTS,@doApplyMPO,{Hl{itN},Hs{itN},Hs{itN+1},Hr{itN+3}});
        
        %using different function to do lanczos diagonalization
        H = reshape(ncon({Hl{itN},Hs{itN},Hs{itN+1},Hr{itN+3}},{[-1,1,-5],[1,-2,2,-6],[2,-3,3,-7],[-4,3,-8]})...
            ,[size(Hl{itN},1)*size(Hs{itN},2)*size(Hs{itN},2)*size(Hr{itN+3},1)...
            ,size(Hl{itN},1)*size(Hs{itN},2)*size(Hs{itN},2)*size(Hr{itN+3},1)]);
        [Aeff, Eeff] = testing_eiglanczos(H,psiGround,[]);
        Eiter(itN,2*itS-1) = Eeff;
        Aeff = reshape(Aeff,[size(Hl{itN},3),size(Hs{itN},4),size(Hs{itN+1},4),size(Hr{itN+3},3)]);
        
        % update A{itN} and B{itN+1} by using Aeff, via SVD
        T = reshape(Aeff,[size(Aeff,1)*size(Aeff,2),size(Aeff,3)*size(Aeff,4)]);
        [U2,S2,V2] = svd(T,'econ');
        N_S_leg = min(min(size(S2)),Nkeep);
        A{itN} = reshape(U2(:,1:N_S_leg),[size(Aeff,1),size(Aeff,2),N_S_leg] );
        sWeight{itN +1} = S2(1:N_S_leg,1:N_S_leg)./sqrt(trace(S2(1:N_S_leg,1:N_S_leg)*S2(1:N_S_leg,1:N_S_leg)'));
        B{itN+1} = reshape(V2(:,1:N_S_leg)',[N_S_leg,size(Aeff,3),size(Aeff,4)] );
        Hr{itN+2} = ncon({Hs{itN+1},Hr{itN+3},B{itN+1},conj(B{itN+1})},{[-2,2,4,1],[3,4,5],[-3,1,5],[-1,2,3]});
    end
        %boundary sigular value
    chil = size(A{1},1);chid = size(A{1},2); chir = size(A{1},3);
    [utemp,stemp,vtemp] = svd(reshape(ncon({A{1},sWeight{2}},{[-1,-2,1],[1,-3]}),[chil,chid*chir]),'econ');
    B{1} = reshape(vtemp',[chil,chid,chir]);
    sWeight{1} = utemp*stemp./sqrt(trace(stemp.^2));
    
    % display informaiton of the sweep
    str = ['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),' (right -> left) : Energy = ',sprintf('%.7g',Eiter(1,2*itS-1))];
    disptime(str);
    % left -> right
    for itN = 1:N-1
        chil = size(B{itN},1); chir = size(B{itN+1},3); chid = size(B{itN},2);
        psiGround = reshape(ncon({sWeight{itN},B{itN},B{itN+1}},{[-1,1],[1,-2,2],[2,-3,-4]}),[chil*chid^2*chir,1]);
        
        %[Aeff, Eeff] = My_eigs_2site(psiGround,OPTS,@doApplyMPO,{Hl{itN},Hs{itN},Hs{itN+1},Hr{itN+3}}); 
        %using different function to do lanczos diagonalization
        H = reshape(ncon({Hl{itN},Hs{itN},Hs{itN+1},Hr{itN+3}},{[-1,1,-5],[1,-2,2,-6],[2,-3,3,-7],[-4,3,-8]})...
            ,[size(Hl{itN},1)*size(Hs{itN},2)*size(Hs{itN},2)*size(Hr{itN+3},1)...
            ,size(Hl{itN},1)*size(Hs{itN},2)*size(Hs{itN},2)*size(Hr{itN+3},1)]);
        [Aeff, Eeff] = testing_eiglanczos(H,psiGround,[]);
        Eiter(itN,2*itS) = Eeff;
        Aeff = reshape(Aeff,[size(Hl{itN},3),size(Hs{itN},4),size(Hs{itN+1},4),size(Hr{itN+3},3)]);
        
        % update A{itN} and B{itN+1} by using Aeff, via SVD
        T = reshape(Aeff,[size(Aeff,1)*size(Aeff,2),size(Aeff,3)*size(Aeff,4)]);
        [U2,S2,V2] = svd(T,'econ');
        N_S_leg = min(min(size(S2)),Nkeep);
        A{itN} = reshape(U2(:,1:N_S_leg),[size(Aeff,1),size(Aeff,2),N_S_leg] );
        sWeight{itN +1} = S2(1:N_S_leg,1:N_S_leg)./sqrt(trace(S2(1:N_S_leg,1:N_S_leg)*S2(1:N_S_leg,1:N_S_leg)'));
        B{itN+1} = reshape(V2(:,1:N_S_leg)',[N_S_leg,size(Aeff,3),size(Aeff,4)] );
        Hl{itN+1} = ncon({Hs{itN},Hl{itN},A{itN},conj(A{itN})},{[4,2,-2,1],[3,4,5],[5,1,-3],[3,2,-1]});
    end
    chil = size(B{itN+1},1); chid = size(B{itN+1},2);chir = size(B{itN+1},3);
    [utemp,stemp,vtemp] = svd(reshape(ncon({B{itN+1},sWeight{itN+1}},{[1,-2,-3],[-1,1]}),[chil*chid,chir,1]),'econ');    
    A{itN+1} = reshape(utemp,[chil,chid,chir]);
    sWeight{itN+2} = (stemp./sqrt(sum(diag(stemp).^2)))*vtemp';
        % display informaiton of the sweep
        str = ['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),' (left -> right) : Energy = ',sprintf('%.7g',Eiter(N-1,2*itS))];
        disptime(str);
end

E0 = Eiter(N,2*itS); % take the last value
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function psi = doApplyMPO(psi,L,M1,M2,R)
        % applies the superblock MPO to the state
        if size(M1)==size(M2)
            psi = reshape(ncon({reshape(psi,[size(L,3),size(M1,4),size(M2,4),size(R,3)]),L,M1,M2,R},...
                {[1,3,5,7],[-1,2,1],[2,-2,4,3],[4,-3,6,5],[-4,6,7]}),[size(L,3)*size(M1,4)*size(M2,4)*size(R,3),1]);
            % when the dmrg sweep reach the boundary the MPO will have three legs
            % that the two MPOs (M1 and M2) will have different size
        elseif size(M1,1) == 1%left boundary Hl == 1
            M1 = permute(M1, [2,3,4,1]);
            psi = reshape(ncon({reshape(psi,[size(M1,3),size(M2,4),size(R,3)]),M1,M2,R},...
                {[1,3,5],[-1,2,1],[2,-2,4,3],[-3,4,5]}),[size(M1,3)*size(M2,4)*size(R,3),1]);
        elseif size(M2,3) == 1%right boundary Hr == 1
            M2 = permute(M2, [1,2,4,3]);
            psi = reshape(ncon({reshape(psi,[size(L,3),size(M1,4),size(M2,3)]),L,M1,M2},...
                {[1,3,5],[-1,2,1],[2,-2,4,3],[4,-3,5]}),[size(L,3)*size(M1,4)*size(M2,3),1]);
        end
    end
    
    
    
    % B = canonForm(M,0); % bring into right-canonical form
% for itN = (N:-1:1)
%     if itN == N
%         % "remove" the right leg (3rd leg) which is dummy by permuting to the last
%         H2 = permute(Hs{itN},[1 2 4 3]);
%         Hr{itN+1} = ncon({B{itN},H2,B{itN}},{[-3,2],[-2,1,2],[-1,1]});
%     else
%         Hr{itN+1} = ncon({Hr{itN+2},B{itN},Hs{itN},B{itN}},{[1,2,3],[-3,5,3],[-2,4,2,5],[-1,4,1]});
%     end
% end