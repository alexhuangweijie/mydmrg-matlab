function [A, B, sWeight, Hl, Hr] = My_boundary(i, itN, A, B, sWeight, Hl, Hs, Hr, Aeff, Nkeep)
% <description>
% this function will using the Aeff to update the A,B,Sweight,Hr(or Hl)
% i is the label of the four condition that will appear in the sweeps
% from right to left: (1) at right boundary (2) at left boundary
% from left to right: (3) at left boundary (4) at right boundary
% written by Huang Weijie (2020.5.25)

if  i == 1
    Aeff = reshape(Aeff,[size(Hl{itN},3),size(Hs{itN},4),size(Hs{itN+1},4)]);
    T = reshape(Aeff,[size(Aeff,1)*size(Aeff,2),size(Aeff,3)]);
    [U2,S2,V2] = svd(T,'econ');
    N_S_leg = min(min(size(S2)),Nkeep);
    A{itN} = reshape(U2(:,1:N_S_leg),[size(Aeff,1),size(Aeff,2),N_S_leg] );
    sWeight{itN +1} = S2(1:N_S_leg,1:N_S_leg)./sqrt(trace(S2(1:N_S_leg,1:N_S_leg)*S2(1:N_S_leg,1:N_S_leg)'));
    B{itN+1} = reshape(V2(:,1:N_S_leg)',[N_S_leg,size(Aeff,3)] );
    Hr{itN+2} = ncon({permute(Hs{itN+1},[1,2,4,3]),B{itN+1},conj(B{itN+1})},{[-2,2,1],[-3,2],[-1,1]});
elseif (i == 2)&&(itN == 1)
    Aeff = reshape(Aeff,[size(Hs{itN},4),size(Hs{itN+1},4),size(Hr{itN+3},3)]);
    T = reshape(Aeff,[size(Aeff,1),size(Aeff,2)*size(Aeff,3)]);
    [U2,S2,V2] = svd(T,'econ');
    N_S_leg = min(min(size(S2)),Nkeep);
    A{itN} = reshape(U2(:,1:N_S_leg),[1,size(Aeff,1),N_S_leg] );
    sWeight{itN +1} = S2(1:N_S_leg,1:N_S_leg)./sqrt(trace(S2(1:N_S_leg,1:N_S_leg)*S2(1:N_S_leg,1:N_S_leg)'));
    B{itN+1} = reshape(V2(:,1:N_S_leg)',[N_S_leg,size(Aeff,2),size(Aeff,3)] );
    Hr{itN+2} = ncon({Hs{itN+1},Hr{itN+3},B{itN+1},conj(B{itN+1})},{[-2,2,4,1],[3,4,5],[-3,1,5],[-1,2,3]});
    %boundary sigular value
    chil = size(A{1},1);chid = size(A{1},2); chir = size(A{1},3);
    [utemp,stemp,vtemp] = svd(reshape(ncon({A{1},sWeight{2}},{[-1,-2,1],[1,-3]}),[chil,chid*chir]),'econ');
    B{1} = reshape(vtemp',[chil,chid,chir]);
    sWeight{1} = utemp*stemp./sqrt(trace(stemp.^2));
    
elseif i == 3
    Aeff = reshape(Aeff,[size(Hs{itN},4),size(Hs{itN+1},4),size(Hr{itN+3},3)]);
    T = reshape(Aeff,[size(Aeff,1),size(Aeff,2)*size(Aeff,3)]);
    [U2,S2,V2] = svd(T,'econ');
    N_S_leg = min(min(size(S2)),Nkeep);
    A{itN} = reshape(U2(:,1:N_S_leg),[size(Aeff,1),N_S_leg] );
    sWeight{itN +1} = S2(1:N_S_leg,1:N_S_leg)./sqrt(trace(S2(1:N_S_leg,1:N_S_leg)*S2(1:N_S_leg,1:N_S_leg)'));
    B{itN+1} = reshape(V2(:,1:N_S_leg)',[N_S_leg,size(Aeff,2),size(Aeff,3)] );
    Hl{itN+1} = ncon({Hs{itN},Hl{itN},A{itN},conj(A{itN})},{[4,2,-2,1],[3,4,5],[5,1,-3],[3,2,-1]});
elseif i == 4
    Aeff = reshape(Aeff,[size(Hl{itN},3),size(Hs{itN},4),size(Hs{itN+1},4)]);
    T = reshape(Aeff,[size(Aeff,1)*size(Aeff,2),size(Aeff,3)]);
    [U2,S2,V2] = svd(T,'econ');
    N_S_leg = min(min(size(S2)),Nkeep);
    A{itN} = reshape(U2(:,1:N_S_leg),[size(Aeff,1),size(Aeff,2),N_S_leg] );
    sWeight{itN +1} = S2(1:N_S_leg,1:N_S_leg)./sqrt(trace(S2(1:N_S_leg,1:N_S_leg)*S2(1:N_S_leg,1:N_S_leg)'));
    B{itN+1} = reshape(V2(:,1:N_S_leg)',[N_S_leg,size(Aeff,3)] );
    Hl{itN+1} = ncon({Hs{itN},Hl{itN},A{itN},conj(A{itN})},{[4,2,-2,1],[3,4,5],[5,1,-3],[3,2,-1]});
    %boundary sigular value
    chil = size(B{itN+1},1); chid = size(B{itN+1},2);chir = size(B{itN+1},3);
    [utemp,stemp,vtemp] = svd(ncon({B{itN+1},sWeight{itN+1}},{[1,-2],[1,-1]}),'econ');
    %A{itN+1} = reshape(utemp,[chil,chid,chir]);
    sWeight{Nsites+1} = (stemp./sqrt(sum(diag(stemp).^2)))*vtemp';
    end
end