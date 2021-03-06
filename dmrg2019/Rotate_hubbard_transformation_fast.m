function [Block1,Block2,T1,T2,Psi_L,Psi_R] = Rotate_hubbard_transformation_fast(A,B,Psi,noftk,transform1,transform2)
Block1 = A;
Block2 = B;
a1 = Block1.basis_size;
b1 = Block2.basis_size;
%左边的BLOCK的处理
Psimatrix = reshape ( Psi , b1 , a1 );
RDM = Psimatrix' * Psimatrix;
%对角化约化密度矩阵 处理左边的块
[ vector1 ,eigenvalues ] = eig ( RDM );
[ ~ , index ] = sort ( diag ( eigenvalues ) ,'descend' );
vector1 = vector1 ( : , index );
numbers_to_keep = min ( size ( vector1 , 1 ) , noftk );
T1 = vector1 ( : , 1 : numbers_to_keep );
Block1.basis_size = numbers_to_keep;
Block1.H  = T1' * Block1.H * T1;
Block1.C = T1' * Block1.C * T1;
Block1.Cdag = T1' * Block1.Cdag * T1;
Block1.N = T1' * Block1.N * T1;
Block1.I = T1' * Block1.I * T1;
%右边的BLOCK的处理
Psimatrix = reshape ( Psi , b1 , a1 );
RDM = Psimatrix * Psimatrix';
[ vector2 , eigenvalues ] = eig( RDM );
[ ~ , index ] = sort( diag ( eigenvalues ) , 'descend' );
vector2 = vector2( : , index );
numbers_to_keep = min( size ( vector2 , 1 ) , noftk );
T2 = vector2( : , 1 : numbers_to_keep );
Block2.basis_size = numbers_to_keep;
Block2.H  = T2' * Block2.H * T2;
Block2.C = T2' * Block2.C * T2;
Block2.Cdag = T2' * Block2.Cdag * T2;
Block2.N = T2' * Block2.N * T2;
Block2.I = T2' * Block2.I * T2;
Psi_L = kron(kron(T1',eye(2)),transform2)*Psi;
Psi_R = kron(kron(transform1,eye(2)),T2')*Psi;
Psi_L = Psi_L./norm(Psi_L);
Psi_R = Psi_R./norm(Psi_R);
end