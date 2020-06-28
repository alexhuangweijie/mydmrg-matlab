%%%%% Ex4.2(a): set B-C link as center of orthogonality
d = 5; % index dimension
A = rand(d,d,d);
B = rand(d,d,d);
C = rand(d,d,d);
Sig = eye(d); % initial link matrix

% generate gauge change matrices
rho1 = ncon({A,A,B,B},{[1,2,3],[1,2,4],[3,5,-1],[4,5,-2]});
rho2 = ncon({C,C},{[-1,1,2],[-2,1,2]});
[u1,d1] = eig(rho1); sq_d1 = sqrt(abs(diag(d1)));
[u2,d2] = eig(rho2); sq_d2 = sqrt(abs(diag(d2)));
X1 = u1*diag(sq_d1)*u1'; X1inv = u1*diag(1./sq_d1)*u1';
X2 = u2*diag(sq_d2)*u2'; X2inv = u2*diag(1./sq_d2)*u2';
% implement gauge change
Bprime = ncon({B,X1inv},{[-1,-2,1],[1,-3]});
Cprime = ncon({X2inv,C},{[-1,1],[1,-2,-3]});
Sig_prime = X1*Sig*X2;
% check result
H0 = ncon({A,B,C},{[-1,-2,1],[1,-3,2],[2,-4,-5]});
H1 = ncon({A,Bprime,Sig_prime,Cprime},{[-1,-2,1],[1,-3,2],[2,3],[3,-4,-5]});
totErr = norm(H0(:) - H1(:)) / norm(H0(:))