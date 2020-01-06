function [ H , H1 , H4 ] = BlockBuilding_hubbard ( H1 , H2 , H3 , H4 ,  J , Bz )
H_Bz = 0;
H1.H  = kron( H1.H , H2.I ) + kron( H1.I , H2.H) + J *( kron( H1.Cdag , H2.C ) +  kron( H1.C , H2.Cdag ) );
H1.C = kron( H1.I , H2.C );
H1.Cdag = kron( H1.I , H2.Cdag );
H1.N =  kron( H1.I , H2.N );
H1.I =  kron( H1.I , H2.I );
H1.basis_size = H1.basis_size * H2.basis_size;
H4.H  =  kron( H3.H , H4.I ) + kron( H3.I , H4.H ) + J *( kron( H3.Cdag , H4.C ) + kron( H3.C , H4.Cdag ) );
H4.C =  kron( H3.C , H4.I );
H4.Cdag = kron( H3.Cdag , H4.I );
H4.N =  kron( H3.N , H4.I );
H4.I =  kron( H3.I , H4.I );
H4.basis_size = H3.basis_size * H4.basis_size;
H = kron( H1.H , H4.I ) + kron( H1.I , H4.H ) + J *( kron( H1.Cdag , H4.C ) +  kron( H1.C , H4.Cdag ) ) + Bz * H_Bz ;
end


