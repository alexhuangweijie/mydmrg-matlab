
SS3 = zeros(N,N);
A = A0;
itp = 180;%set 90 degree
% % rotation
U_1_rotation_bond = zeros(2,2,2,2);
%rotation angle of the U(1) rotation
U_1_rotation_bond(1,1,1,1) = 1;
U_1_rotation_bond(1,2,1,2) = exp(-2*pi/360*itp*1i);
U_1_rotation_bond(2,1,2,1) = exp(-2*pi/360*itp*1i);
U_1_rotation_bond(2,2,2,2) = 1;
U_1_rotation_bond = permute(U_1_rotation_bond,[1,4,3,2]);
U_1 = U_1_rotation_bond;
%G(itN_1, itN_2)

% calculate all the Gij pairs 
for itN_1 = 1:N-1
    if itN_1 == 1
        %Sx@site1
        tensors_N1 = {A{1},A{1},U_1,conj(U_1),Z,Sx*2};
        connects_N1 = {[6,3,-6],[6,5,-1],[2,-4,-5,3],[4,-3,-2,5],[1,2],[4,1]};
        cont_order_N1 = [1, 2, 5, 3, 6, 4];
        T1 = ncon(tensors_N1, connects_N1 ,cont_order_N1);
    else
        %empty @site1
        tensors_N2 = {A{1},U_1,conj(U_1),A{1}};
        connects_N2 = {[1,2,-6],[3,-4,-5,2],[3,-3,-2,4],[1,4,-1]};
        cont_order_N2 = [4, 2, 1, 3];
        T1 = ncon(tensors_N2, connects_N2 ,cont_order_N2);
        if itN_1 > 2
            for itN_3 = 2 : itN_1-1
                %filling in the blank
                tensors_N1 = {T1,U_1,conj(U_1),A{itN_3},A{itN_3}};
                connects_N1 = {[7,1,2,4,3,6],[5,-4,-5,4],[5,-3,-2,2],[6,3,-6],[7,1,-1]};
                cont_order_N1 = [7, 1, 5, 3, 6, 2, 4];
                T1 = ncon(tensors_N1, connects_N1 ,cont_order_N1);
            end
        end
        %Sx@site itN_1
        tensors_N3 = {T1,A{itN_1},U_1,conj(U_1),A{itN_1},Z,Sx*2};
        connects_N3 = {[1,4,6,5,3,2],[2,3,-6],[9,-4,-5,5],[7,-3,-2,6],[1,4,-1],[8,9],[7,8]};
        cont_order_N3 = [8, 7, 1, 4, 3, 2, 5, 6, 9];
        T1 = ncon(tensors_N3, connects_N3 ,cont_order_N3);
    end
    for itN_2 = itN_1+1 : N
        T = T1;
        if itN_2 > itN_1+1
            for itN_4 = itN_1+1:itN_2-1
                %using Z operator to filling the blank
                tensors_N2 = {T,A{itN_4},U_1,conj(U_1),A{itN_4},Z,};
                connects_N2 = {[1,5,6,4,3,2],[2,3,-6],[8,-4,-5,4],[7,-3,-2,6],[1,5,-1],[7,8]};
                cont_order_N2 = [7, 1, 5, 3, 2, 6, 4, 8];
                T = ncon(tensors_N2, connects_N2 ,cont_order_N2);
            end
        end
        if itN_2 == N
            %Sy@site end
            tensors_N4 = {Sy*2,T,A{N},A{N}};
            connects_N4 = {[6,5],[2,3,6,5,4,1],[1,4],[2,3]};
            cont_order_N4 = [2, 3, 6, 5, 4, 1];
            T = ncon(tensors_N4, connects_N4 ,cont_order_N4);
        else
            %Sy@site itN_2
            tensors_N4 = {T,A{itN_2},U_1,Sy*2,conj(U_1),A{itN_2}};
            connects_N4 = {[1,4,5,6,3,2],[2,3,-6],[7,-4,-5,6],[8,7],[8,-3,-2,5],[1,4,-1]};
            cont_order_N4 = [8, 7, 1, 4, 3, 2, 5, 6];
            T = ncon(tensors_N4, connects_N4 ,cont_order_N4);
            if itN_2 < N-1
                for itN_5 = itN_2+1 : N-1
                    %filling the blank with empty
                    tensors_N1 = {T,U_1,conj(U_1),A{itN_5},A{itN_5}};
                    connects_N1 = {[7,1,2,4,3,6],[5,-4,-5,4],[5,-3,-2,2],[6,3,-6],[7,1,-1]};
                    cont_order_N1 = [7, 1, 5, 3, 6, 2, 4];
                    T = ncon(tensors_N1, connects_N1 ,cont_order_N1);
                end
            end
            %empty@the end
            tensors_N3 = {T,A{N},A{N}};
            connects_N3 = {[1,2,5,5,4,3],[3,4],[1,2]};
            cont_order_N3 = [5, 1, 2, 4, 3];
            T = ncon(tensors_N3, connects_N3 ,cont_order_N3);
        end
        SS3(itN_1,itN_2) = T;
    end
end

figure(2)
[X, Y] = meshgrid(1:N, 1:N);
meshz(X, Y, abs(SS3))
zlabel('$|G_{ij}|$','Interpreter','latex','fontsize',20);
xlabel('j');
ylabel('i');







