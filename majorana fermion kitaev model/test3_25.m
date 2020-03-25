clear
j.t=1;
j.delta=0;
j.mu=0;
j.u=0.5;
Bz=0;

%Psi = [1;5;3;14];
BLOCK1.basis_size = 2;
BLOCK1.I = eye(2);
BLOCK1.C = [0,0;1,0];%lambda1
BLOCK1.Cdag = [0,1;0,0];%lambda2
BLOCK1.H = [0,0;0,0];
BLOCK1.N = [0,0;0,1];
BLOCK2.basis_size = 2;
BLOCK2.I = eye(2);
BLOCK2.C = [0,0;1,0];%lambda1
BLOCK2.Cdag = [0,1;0,0];%lambda2
BLOCK2.H = [0,0;0,0];
BLOCK2.N = [0,0;0,1];

Block_L = BLOCK1;
Block_R = BLOCK1;

% for i = 1:2
%     [Block_L,Block_R,H_super]=BlockBuilding_hubbard(Block_L,BLOCK2,BLOCK2,Block_R,j,Bz);
% end
D = 2;
BLOCK2.basis_size = D;
BLOCK2.I = eye(D);
BLOCK2.C = rand(D,D);%lambda1
BLOCK2.Cdag = rand(D,D);%lambda2
BLOCK2.H = zeros(D,D);
BLOCK2.N = ones(D,D);
Block_L = BLOCK2;
Block_R = BLOCK2;
[Block_L,Block_R,H_super]=BlockBuilding_hubbard(Block_L,BLOCK2,BLOCK2,Block_R,j,Bz);
[PsiL,T1,EnergyL] = svd_L(Block_L ,Block_R , j , 0);
[PsiR,T2,EnergyR] = svd_R(Block_L ,Block_R , j , 0);
%H1 = H_super; H1(1,1) = 10;H_super = H1;
H = H_super'*H_super;
[a,b]=eig(H);
c = diag(b);
c(1:3)
%det(H)
%H(1,1) = 100000;
det(H)
