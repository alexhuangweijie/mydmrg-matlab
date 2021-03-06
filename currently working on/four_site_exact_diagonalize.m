



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%testing the U_1 rotation method by using four site hamiltonian with exact diagonalize.
%the phase of each state will be exp(-i*a*n), a is the phase of the bond
%between state |0> and |1>, n is the number of the bond between |0> and |1>
% hwj/2020/9/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comparing the way adding the U_1 by applying many operator to the wave
%function with the exact operator.
% hwj/2020/9/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa = [];c = [];
% for ii = 1:11
U = 0.5;
mu = 0.5;
t = 1;
I = eye(2);
Sx = 0.5*[0,1;1,0];
Sz = 0.5*[1,0;0,-1];
Sy = 0.5*[0,1;-1,0];
Z = [-1,0;0,1];
%finding ground state
H = mu * ( kron(Sz,eye(32)) + kron(kron(eye(2),Sz),eye(16)) + kron(kron(eye(4),Sz),eye(8)) + kron(kron(eye(8),Sz),eye(4)) + kron(kron(eye(16),Sz),eye(2)) + kron(eye(32),Sz) )...
    + 4*U * ( kron(kron(Sz,Sz),eye(16)) + kron(kron(kron(eye(2),Sz),Sz),eye(8)) + kron(kron(kron(eye(4),Sz),Sz),eye(4)) + kron(kron(kron(eye(8),Sz),Sz),eye(2)) + kron(kron(eye(16),Sz),Sz) )...
    - 4*t * ( kron(kron(Sx,Sx),eye(16)) + kron(kron(kron(eye(2),Sx),Sx),eye(8)) + kron(kron(kron(eye(4),Sx),Sx),eye(4)) + kron(kron(kron(eye(8),Sx),Sx),eye(2)) + kron(kron(eye(16),Sx),Sx) );
[a, b] = eigs(H,1,'sr');
c = [];
for i = 1:64
    if abs(a(i))>0.001
        c = [c;i];
    end
end
res = dec2bin(64-c);
aa = dec2bin(64-linspace(1,64,64));
%adding phase and calculate Gij
U_1_angle = [
    0
    1
    2
    1
    2
    3
    2
    1
    2
    3
    4
    3
    2
    3
    2
    1
    2
    3
    4
    3
    4
    5
    4
    3
    2
    3
    4
    3
    2
    3
    2
    1
    1
    2
    3
    2
    3
    4
    3
    2
    3
    4
    5
    4
    3
    4
    3
    2
    1
    2
    3
    2
    3
    4
    3
    2
    1
    2
    3
    2
    1
    2
    1
    0
    ];

index = c;phase = zeros(360,3);
for angle = 1:360
    U_1 = eye(64,64);
    for i = 1:64
        U_1(i,i) = exp(-1i/360*pi*angle*10*U_1_angle(i));
    end
    MPO = cell(6,6);
    MPO{1,2} = kron(kron(Sx*2*Z,Sy*2),eye(16));
    MPO{1,3} = kron(kron(kron(Sx*2*Z,Z),Sy*2),eye(8));
    MPO{1,4} = kron(kron(kron(kron(Sx*2*Z,Z),Z),Sy*2),eye(4));
    MPO{1,5} = kron(kron(kron(kron(kron(Sx*2*Z,Z),Z),Z),Sy*2),eye(2));
    MPO{1,6} = kron(kron(kron(kron(kron(Sx*2*Z,Z),Z),Z),Z),Sy*2);
    MPO{2,3} = kron(kron(kron(eye(2),Sx*2*Z),Sy*2),eye(8));
    MPO{2,4} = kron(kron(kron(kron(eye(2),Sx*2*Z),Z),Sy*2),eye(4));
    MPO{2,5} = kron(kron(kron(kron(kron(eye(2),Sx*2*Z),Z),Z),Sy*2),eye(2));
    MPO{2,6} = kron(kron(kron(kron(kron(eye(2),Sx*2*Z),Z),Z),Z),Sy*2);
    MPO{3,4} = kron(kron(kron(eye(4),Sx*2*Z),Sy*2),eye(4));
    MPO{3,5} = kron(kron(kron(kron(eye(4),Sx*2*Z),Z),Sy*2),eye(2));
    MPO{3,6} = kron(kron(kron(kron(eye(4),Sx*2*Z),Z),Z),Sy*2);
    MPO{4,5} = kron(kron(kron(eye(8),Sx*2*Z),Sy*2),eye(2));
    MPO{4,6} = kron(kron(kron(eye(8),Sx*2*Z),Z),Sy*2);
    MPO{5,6} = kron(kron(eye(16),Sx*2*Z),Sy*2);
    MPO{1,1} = kron(Sz*2,eye(32));
    MPO{2,2} = kron(kron(eye(2),Sz*2),eye(16));
    MPO{3,3} = kron(kron(eye(4),Sz*2),eye(8));
    MPO{4,4} = kron(kron(eye(8),Sz*2),eye(4));
    MPO{5,5} = kron(kron(eye(16),Sz*2),eye(2));
    MPO{6,6} = kron(eye(32),Sz*2);
    
    G = zeros(6,6);GG = zeros(6,6);
    for index_i = 1:6
        for index_j = index_i:6
            if MPO{index_i,index_j} == 0
                MPO{index_i,index_j} = eye(64,64);
            end
            G(index_i, index_j) = a'*U_1'*MPO{index_i,index_j}*U_1*a;
            GG = abs(G);
        end
    end
    
    
    GG = abs(G);
    phase(angle,1) = G(1,6);
    phase(angle,2) = G(1,2);
    phase(angle,3) = G(2,3);
end
plot(real(phase))
% [X, Y] = meshgrid(1:6, 1:6);
% meshz(X, Y, GG)
% zlabel('$|G_{ij}|$','Interpreter','latex','fontsize',20);
% xlabel('j');
% ylabel('i');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%operator testing in 6 sites 
itp = 5;
U_1_rotation_bond = zeros(2,2,2,2);
%rotation angle of the U(1) rotation
U_1_rotation_bond(1,1,1,1) = 1;
U_1_rotation_bond(1,2,1,2) = exp(1i*pi/360*itp);
U_1_rotation_bond(2,1,2,1) = exp(1i*pi/360*itp);
U_1_rotation_bond(2,2,2,2) = 1;

tensors = {U_1_rotation_bond, U_1_rotation_bond, U_1_rotation_bond,...
           U_1_rotation_bond, U_1_rotation_bond};
connects = {[-1,1,-11,-12],[-2,2,-10,1],[-3,3,-9,2],[-4,4,-8,3],[-5,-6,-7,4]};
a = ncon(tensors, connects);
%a = reshape(a,[64,64]);
%exact result
    for i = 1:64
        U_1(i,i) = exp(1i/360*pi*itp*U_1_angle(i));
    end
C = find(a);
D = find(U_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%testing operator in 2 sites
itp = 5;
U_1_rotation_bond = zeros(2,2,2,2);
%rotation angle of the U(1) rotation
U_1_rotation_bond(1,1,1,1) = 1;
U_1_rotation_bond(1,2,1,2) = exp(1i*pi/360*itp);
U_1_rotation_bond(2,1,2,1) = exp(1i*2*pi/360*itp);
U_1_rotation_bond(2,2,2,2) = 1;
tensors = {U_1_rotation_bond,  U_1_rotation_bond};
connects = {[-1,1,-5,-4],[-2,-3,-6,1]};
a = ncon(tensors, connects);












