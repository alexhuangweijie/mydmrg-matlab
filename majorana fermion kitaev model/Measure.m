function [ m1, m2] = Measure ( H1 , H2 ,  H3 , H4 , PsiL , PsiR )
if (size(PsiL,1)==H1.basis_size*H4.basis_size)
    I1 = eye( H1.basis_size / H2.basis_size );
    I2 = eye( H4.basis_size / H3.basis_size );
else
    I1 = eye( H1.basis_size );
    I2 = eye( H4.basis_size );
end 
    %m1 =  Psi' * kron( kron( I1,kron(H2.N, H3.I)), I2) * Psi;
    m1 =  PsiL' * mykron(kron( I1,H2.N), kron( H3.I,I2), PsiR);
    %m2 =  Psi' * kron( kron( I1,kron(H2.I, H3.N)), I2) * Psi;
    m2 =  PsiL' * mykron( kron( I1,H2.I), kron(H3.N,I2), PsiR);
    m1 = (m1*1i+1)*0.5;
    m2 = (m2*1i+1)*0.5;    
end

    function [e] = mykron(a,b,p)
        d = reshape(p,length(b),length(a));
        e = b*d*(a.');
        e = e(:);
        %e2 = e.';
        %e3 = e2(:);
    end