clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  testing hamitonian from miaojianjian PRL using Sz Sx %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding the U(1) rotation to every bonds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the parameters are set in this page
%system parameter
format long
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
Hs{end} = Hloc1(:,:,1,:);%rotation_bond
%Hs{1} = Hloc1(end,:,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A0,E0,Eiter0] = DMRG_2site (Hs,Nkeep,Nsweep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding the U(1) rotation to every bonds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcultate G1L for different angle of phase 
% A = A0;
% phase = zeros(360,3);
% for itp = 1:360%calcultate G1L for different angle of phase 
%     % % rotation
%     U_1_rotation_bond = zeros(2,2,2,2);
%     %rotation angle of the U(1) rotation
%     U_1_rotation_bond(1,1,1,1) = 1;
%     U_1_rotation_bond(1,2,1,2) = exp(-2*pi/360*itp*1i);
%     U_1_rotation_bond(2,1,2,1) = exp(-2*pi/360*itp*1i);
%     U_1_rotation_bond(2,2,2,2) = 1;
%     U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
% 
%     tensors_N1 = {A{1},A{2},A{1},A{2},U_1_rotation_bond,conj(U_1_rotation_bond),Z,Sz*2,};
%     connects_N1 = {[7,1,5],[5,2,-4],[7,8,6],[6,9,-1],[4,-3,2,1],[10,-2,9,8],[3,4],[10,3]};
%     cont_order_N1 = [3, 4, 6, 8, 9, 5, 1, 2, 7, 10];
%     T4 = ncon(tensors_N1, connects_N1 ,cont_order_N1);
%     for itp4 = 3:N-1
%         tensors_N2 = {T4,A{itp4},A{itp4},Z,U_1_rotation_bond,conj(U_1_rotation_bond),};
%         connects_N2 = {[6,7,4,5],[6,8,-1],[5,3,-4],[1,2],[2,-3,3,4],[1,-2,8,7]};
%         cont_order_N2 = [2, 6, 7, 8, 4, 1, 5, 3];
%         T4 = ncon(tensors_N2, connects_N2, cont_order_N2);
%     end
%     tensors_N3 = {T4,U_1_rotation_bond,Z,Sy*2,A{N},A{N},conj(U_1_rotation_bond)};
%     connects_N3 = {[6,9,4,5],[1,2,3,4],[7,1],[8,2],[5,3],[6,10],[7,8,10,9]};
%     cont_order_N3 = [8, 7, 1, 2, 3, 10, 6, 9, 4, 5];
%     T4 = ncon(tensors_N3,connects_N3,cont_order_N3);
%     phase(itp, 1) = T4;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %calculate G(myi,myi + 1)
% my_i  = 30;
% A = A0;
% itp = 90;
% for index_myi = 1:2
%     my_i = 30 + index_myi;
%     %for itp = 1:360%calcultate G1L for different angle of phase
%         % % rotation
%         U_1_rotation_bond = zeros(2,2,2,2);
%         %rotation angle of the U(1) rotation
%         U_1_rotation_bond(1,1,1,1) = 1;
%         U_1_rotation_bond(1,2,1,2) = exp(-2*pi/360*itp*1i);
%         U_1_rotation_bond(2,1,2,1) = exp(-2*pi/360*itp*1i);
%         U_1_rotation_bond(2,2,2,2) = 1;
%         U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
%         U_1 = U_1_rotation_bond;
%         
%         tensors_N1 = {A{2},U_1,conj(U_1),A{2},A{1},A{1}};
%         connects_N1 = {[4,1,-4],[8,-3,1,6],[8,-2,2,5],[3,2,-1],[7,6,4],[7,5,3]};
%         cont_order_N1 = [5, 2, 3, 6, 4, 1, 8, 7];
%         T = ncon(tensors_N1, connects_N1 ,cont_order_N1);%first with out operator
%         for itN1 = 3:my_i-1
%             tensors_N2 = {U_1,conj(U_1),T,A{itN1},A{itN1}};
%             connects_N2 = {[7,-3,4,1],[7,-2,3,2],[5,2,1,6],[6,4,-4],[5,3,-1]};
%             cont_order_N2 = [6, 4, 1, 7, 2, 5, 3];
%             T = ncon(tensors_N2, connects_N2 ,cont_order_N2);%filing the blank with empty operator
%         end
%         tensors_N3 = {A{my_i},U_1,conj(U_1),T,A{my_i},Z,Sx*2};
%         connects_N3 = {[8,1,-4],[2,-3,1,6],[4,-2,7,5],[9,5,6,8],[9,7,-1],[3,2],[4,3]};
%         cont_order_N3 = [3, 2, 8, 1, 6, 5, 4, 9, 7];
%         T = ncon(tensors_N3, connects_N3 ,cont_order_N3);%my_i
%         tensors_N4 = {U_1,conj(U_1),T,A{my_i+1},A{my_i+1},Sy*2,};
%         connects_N4 = {[8,-3,4,1],[5,-2,3,2],[6,2,1,7],[7,4,-4],[6,3,-1],[5,8]};
%         cont_order_N4 = [6, 8, 3, 2, 1, 5, 4, 7];
%         T = ncon(tensors_N4, connects_N4 ,cont_order_N4);%my_i+1
%         for itN2 = my_i +2 : N-1
%             tensors_N5 = {U_1,conj(U_1),T,A{itN2},A{itN2}};
%             connects_N5 = {[7,-3,4,1],[7,-2,3,2],[5,2,1,6],[6,4,-4],[5,3,-1]};
%             cont_order_N5 = [6, 4, 1, 7, 2, 5, 3];
%             T = ncon(tensors_N5, connects_N5 ,cont_order_N5);%to the end
%         end
%         tensors_N5 = {U_1,conj(U_1),T,A{N},A{N},};
%         connects_N5 = {[7,8,4,1],[7,8,3,2],[5,2,1,6],[6,4,-2],[5,3,-1]};
%         cont_order_N5 = [6, 4, 1, 7, 8, 2, 5, 3];
%         T = ncon(tensors_N5, connects_N5 ,cont_order_N5);%the end
%         phase(itp, index_myi + 1) = T;
%     %end
% end
