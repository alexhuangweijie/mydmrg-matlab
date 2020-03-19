function [ m1, m2] = Measure_hubbard_time ( H1 , H2 ,  H3 , H4 , Psi )
I1 = eye( H1.basis_size );
I2 = eye( H4.basis_size );
m1 =  Psi' * kron( kron( I1,kron(H2.N, H3.I)), I2) * Psi;
m2 =  Psi' * kron( kron( I1,kron(H2.I, H3.N)), I2) * Psi;
end

