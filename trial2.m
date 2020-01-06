
function [A] = trial2(H1,H2,H3,H4,J,C)
% J = 1;
% BLOCK.basis_size = 2;
% BLOCK.I = eye(2);
% BLOCK.C = [0,1;0,0];
% BLOCK.Cdag = [0,0;1,0];
% BLOCK.H = [0,0;0,0];
% BLOCK.N = [0,0;0,1];
% H1 = BLOCK; H2 = BLOCK; H3 = BLOCK; H4 = BLOCK;
S1 = H1.basis_size;
S2 = H2.basis_size;
S3 = H3.basis_size;
S4 = H4.basis_size;
% C = rand(2*2*2*2,1);

%分成七次算 H1 H2 H3 H4 c1c2 c2c3 c3c4
%第一次H1
  H = H1.H;
  m1 = 1;
  m2 = S1;
  m3 = S2*S3*S4;
  A1 = zeros(m2*m3,1);
for index_1_1 = 1:m1
    for index_1_2 = 1:m2
        for  index_1_3 = 1:m3
            for j_index_1 = 1:m2
                A1((index_1_1-1)*m2*m3 + (index_1_2-1)*m3 + index_1_3,1) = A1((index_1_1-1)*m2*m3 + (index_1_2-1)*m3 + index_1_3,1)+H(index_1_2,j_index_1)*C((index_1_1-1)*m2*m3 + (j_index_1-1)*m3 + index_1_3,1);
            end
        end
    end
end
%第二次 H2
  H = H2.H;
  m1 = S1;
  m2 = S2;
  m3 = S3*S4;
  A2 = zeros(m1*m2*m3,1);
for index_2_1 = 1:m1
    for index_2_2 = 1:m2
        for  index_2_3 = 1:m3
            for j_index_2 = 1:m2
                A2((index_2_1-1)*m2*m3 + (index_2_2-1)*m3 + index_2_3,1) = A2((index_2_1-1)*m2*m3 + (index_2_2-1)*m3 + index_2_3,1)+H(index_2_2,j_index_2)*C((index_2_1-1)*m2*m3 + (j_index_2-1)*m3 + index_2_3,1);
            end
        end
    end
end
%第三次H3
  H = H3.H;
  m1 = S1*S2;
  m2 = S3;
  m3 = S4;
  A3 = zeros(m1*m2*m3,1);
for index_3_1 = 1:m1
    for index_3_2 = 1:m2
        for  index_3_3 = 1:m3
            for j_index_3 = 1:m2
                A3((index_3_1-1)*m2*m3 + (index_3_2-1)*m3 + index_3_3,1) = A3((index_3_1-1)*m2*m3 + (index_3_2-1)*m3 + index_3_3,1)+H(index_3_2,j_index_3)*C((index_3_1-1)*m2*m3 + (j_index_3-1)*m3 + index_3_3,1);
            end
        end
    end
end

%第四次H4
 
  H = H4.H;
  m1 = S1*S2*S3;
  m2 = S4;
  m3 = 1;
  A4 = zeros(m1*m2,1);
for index_4_1 = 1:m1
    for index_4_2 = 1:m2
        for  index_4_3 = 1:m3
            for j_index_4 = 1:m2
                A4((index_4_1-1)*m2*m3 + (index_4_2-1)*m3 + index_4_3,1) = A4((index_4_1-1)*m2*m3 + (index_4_2-1)*m3 + index_4_3,1)+H(index_4_2,j_index_4)*C((index_4_1-1)*m2*m3 + (j_index_4-1)*m3 + index_4_3,1);
            end
        end
    end
end

%第五次c1c2
  H = J *( kron( H1.Cdag , H2.C ) +  kron( H1.C , H2.Cdag ) );
  m1 = 1;
  m2 = S1*S2;
  m3 = S3*S4;
  A5 = zeros(m2*m3,1);
for index_5_1 = 1:m1
    for index_5_2 = 1:m2
        for  index_5_3 = 1:m3
            for j_index_5 = 1:m2
                A5((index_5_1-1)*m2*m3 + (index_5_2-1)*m3 + index_5_3,1) = A5((index_5_1-1)*m2*m3 + (index_5_2-1)*m3 + index_5_3,1)+H(index_5_2,j_index_5)*C((index_5_1-1)*m2*m3 + (j_index_5-1)*m3 + index_5_3,1);
            end
        end
    end
end
%第六次c2c3
  H = J *( kron( H2.Cdag , H3.C ) +  kron( H2.C , H3.Cdag ) );
  m1 = S1;
  m2 = S2*S3;
  m3 = S4;
  A6 = zeros(m1*m2*m3,1);
for index_6_1 = 1:m1
    for index_6_2 = 1:m2
        for  index_6_3 = 1:m3
            for j_index_6 = 1:m2
                A6((index_6_1-1)*m2*m3 + (index_6_2-1)*m3 + index_6_3,1) = A6((index_6_1-1)*m2*m3 + (index_6_2-1)*m3 + index_6_3,1)+H(index_6_2,j_index_6)*C((index_6_1-1)*m2*m3 + (j_index_6-1)*m3 + index_6_3,1);
            end
        end
    end
end
%第七次c3c4
  H = J *( kron( H3.Cdag , H4.C ) +  kron( H3.C , H4.Cdag ) );
  m1 = S1*S2;
  m2 = S3*S4;
  m3 = 1;
  A7 = zeros(m1*m2*m3,1);
for index_1_1 = 1:m1
    for index_1_2 = 1:m2
        for  index_1_3 = 1:m3
            for j_index_1 = 1:m2
                A7((index_1_1-1)*m2*m3 + (index_1_2-1)*m3 + index_1_3,1) = A7((index_1_1-1)*m2*m3 + (index_1_2-1)*m3 + index_1_3,1)+H(index_1_2,j_index_1)*C((index_1_1-1)*m2*m3 + (j_index_1-1)*m3 + index_1_3,1);
            end
        end
    end
end

A = A1+A2+A3+A4+A5+A6+A7;
% [a,~,~] = BlockBuilding_hubbard(H1,H2,H3,H4,J,0);
% a1 = a*C;
end