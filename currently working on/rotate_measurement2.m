function T = rotate_measurement2(A, N, U_1, itN_1, itN_2)
% this function will doing the U(1) rotation to the wave function getting
% from the variational dmrg algorithm, and then calculate the correlation
% function
% 
% 
% input : 
% A : wave function in mps language, which is a cell, A{i} will be the tensor 
% of site i
% 
% N : this is the number of the sites in the chain
% 
% itN_1 and itN_2 : G(itN_1, itN_2) will be calculated , which is Sx*2 in
% itN_1, Sy*2 in site itN_2.
% the Z will be used as thess are fermions which are antisysmetric.
% 
% output : T = G(itN_1, itN_2)  
% Updated by W.huang( September 10, 2020): the Z operator should be added
% to Sy instead of Sx.
Z = [-1,0;0,1];%Z
Sx = 0.5*[0,1;1,0];%Sx
Sz = 0.5*[1,0;0,-1];%Sz
Sy = 0.5*[0,1;-1,0];% i*Sy

if itN_1 == 1
    %Sx@site1
    tensors_N1 = {A{1},A{1},U_1,conj(U_1),Sx*2};
    connects_N1 = {[4,1,-6],[4,3,-1],[5,-4,1,-5],[2,-3,3,-2],[2,5]};
    cont_order_N1 = [5, 3, 1, 4, 2];
    T1 = ncon(tensors_N1, connects_N1 ,cont_order_N1);
else
    %empty @site1
    tensors_N2 = {A{1},U_1,conj(U_1),A{1}};
    connects_N2 = {[1,2,-6],[3,-4,2,-5],[3,-3,4,-2],[1,4,-1]}; 
    cont_order_N2 = [4, 2, 1, 3];
    T1 = ncon(tensors_N2, connects_N2 ,cont_order_N2);
    if itN_1 > 2
        for itN_3 = 2 : itN_1-1
            %filling in the blank
            tensors_N1 = {T1,U_1,conj(U_1),A{itN_3},A{itN_3}};
            connects_N1 = {[6,1,7,3,2,5],[4,-4,3,-5],[4,-3,7,-2],[5,2,-6],[6,1,-1]};
            cont_order_N1 = [6, 1, 4, 2, 5, 7, 3]; 
            T1 = ncon(tensors_N1, connects_N1 ,cont_order_N1);
        end
    end
    %Sx@site itN_1
    tensors_N3 = {T1,A{itN_1},U_1,conj(U_1),A{itN_1},Sx*2};
    connects_N3 = {[1,4,6,5,3,2],[2,3,-6],[8,-4,5,-5],[7,-3,6,-2],[1,4,-1],[7,8]}; 
    cont_order_N3 = [7, 1, 4, 3, 2, 6, 5, 8];
    T1 = ncon(tensors_N3, connects_N3 ,cont_order_N3);
end
T = T1;
if itN_2 > itN_1+1
    for itN_4 = itN_1+1:itN_2-1
        %using Z operator to filling the blank
        tensors_N2 = {T,A{itN_4},U_1,conj(U_1),A{itN_4},Z};
        connects_N2 = {[1,5,6,4,3,2],[2,3,-6],[8,-4,4,-5],[7,-3,6,-2],[1,5,-1],[7,8]};
        cont_order_N2 = [7, 1, 5, 3, 2, 6, 4, 8];
        T = ncon(tensors_N2, connects_N2 ,cont_order_N2);
    end
end
if itN_2 == N
    %Sy@site end
    tensors_N4 = {Sy*2,T,A{N},A{N},Z};
    connects_N4 = {[7,5],[2,3,6,5,4,1],[1,4],[2,3],[6,7]};
    cont_order_N4 = [4, 1, 2, 3, 5, 7, 6];
    T = ncon(tensors_N4, connects_N4 ,cont_order_N4);
else
    %Sy@site itN_2
    tensors_N4 = {T,A{itN_2},U_1,Sy*2,conj(U_1),A{itN_2},Z}; 
    connects_N4 = {[1,4,5,6,3,2],[2,3,-6],[7,-4,6,-5],[9,7],[8,-3,5,-2],[1,4,-1],[8,9]};
    cont_order_N4 = [9, 8, 1, 4, 3, 2, 6, 5, 7];
    T = ncon(tensors_N4, connects_N4 ,cont_order_N4);
    if itN_2 < N-1
        for itN_5 = itN_2+1 : N-1
            %filling the blank with empty
            tensors_N1 = {T,U_1,conj(U_1),A{itN_5},A{itN_5}};
            connects_N1 = {[6,1,7,3,2,5],[4,-4,3,-5],[4,-3,7,-2],[5,2,-6],[6,1,-1]}; 
            cont_order_N1 = [6, 1, 4, 2, 5, 7, 3]; 
            T = ncon(tensors_N1, connects_N1 ,cont_order_N1);
        end
    end
    %empty@the end
    tensors_N3 = {T,A{N},A{N}};
    connects_N3 = {[1,2,5,5,4,3],[3,4],[1,2]};
    cont_order_N3 = [5, 1, 2, 4, 3];
    T = ncon(tensors_N3, connects_N3 ,cont_order_N3);
end
%         SS3(itN_1,itN_2) = T;
end









