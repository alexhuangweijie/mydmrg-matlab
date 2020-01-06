function [ Block1 , Block2 ] = Rotate_hubbard ( A , B , Psi , noftk   )
Block1 = A;
Block2 = B;
a1 = Block1.basis_size;
b1 = Block2.basis_size;
%左边的BLOCK的处理
Psimatrix = reshape ( Psi , b1 , a1 );
RDM = Psimatrix' * Psimatrix;
%对角化约化密度矩阵
[ vector ,eigenvalues ] = eig ( RDM );
[ ~ , index ] = sort ( diag ( eigenvalues ) ,'descend' );
vector = vector ( : , index );
numbers_to_keep = min ( size ( vector , 1 ) , noftk );
T = vector ( : , 1 : numbers_to_keep );
Block1.basis_size = numbers_to_keep;
Block1.H  = T' * Block1.H * T;
Block1.C = T' * Block1.C * T;
Block1.Cdag = T' * Block1.Cdag * T;
Block1.N = T' * Block1.N * T;
Block1.I = T' * Block1.I * T;
%右边的BLOCK的处理
Psimatrix = reshape ( Psi , b1 , a1 );
RDM = Psimatrix * Psimatrix';
[ vector , eigenvalues ] = eig( RDM );
[ ~ , index ] = sort( diag ( eigenvalues ) , 'descend' );
vector = vector( : , index );
numbers_to_keep = min( size ( vector , 1 ) , noftk );
T = vector( : , 1 : numbers_to_keep );
Block2.basis_size = numbers_to_keep;
Block2.H  = T' * Block2.H * T;
Block2.C = T' * Block2.C * T;
Block2.Cdag = T' * Block2.Cdag * T;
Block2.N = T' * Block2.N * T;
Block2.I = T' * Block2.I * T;
end