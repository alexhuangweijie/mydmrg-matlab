function [V,E] = testing_eiglanczos (H,Vinit,Vorth)
% default parameter
N = 20;
minH = 1e-10;
Vinit = Vinit(:);
 if isempty(Vorth), Vorth = zeros(size(H,1),0); end
H = (H+H')/2; % make it Hermitian (to rule out possible numerical noise)
U = [Vorth,Vinit];
N0 = size(U,2);
N = min(N,size(H,2)-N0);
% basis vectors; by adding 0's, we assign the memory space for the vectors to appear later
U = [U,zeros(size(U,1),N)];
ff = zeros(N,1); % 1st diagonal
gg = zeros(N+1,1); % main diagonal

for it1 = (1:(N+1))
    v = H*U(:,N0+it1-1);
    gg(it1) = U(:,N0+it1-1)'*v; % diagonal element; "on-site energy"
    
    if it1 < (N+1)
        % orthonormalize v against U
        Utmp = U(:,(1:(N0+it1-1))); % basis vectors up to the last iteration
        v = v - Utmp*(Utmp'*v); % project out the component which is parallel to the column vectors of Utmp
        v = v - Utmp*(Utmp'*v); % do twice, to further suppress numerical noise
        
        nv = norm(v); % norm of the vector
        
        if nv > minH
            ff(it1) = nv; % 1st diagonal element; "hopping amplitude"
            U(:,N0+it1) = v/ff(it1);
        else
            % stop iteration; truncate ff, gg, and U
            ff(it1:end) = [];
            gg(it1+1:end) = [];
            U(:,(N0+it1):end) = [];
            break;
        end
    end
end

Heff = diag(ff,1);
Heff = Heff + Heff' + diag(gg);

% Use eigs, not eig, since we need only one eigenvector with the lowest eigenvalue
[V,~] = eigs(Heff,1,'sa'); % 'sa': smallest real
V = U(:,(end-size(V,1)+1):end)*V; % basis transformation

% the eigenvalue of Heff is inaccurate measure of the ground-state E, since
% Heff is the matrix constrained to the Krylov subspace.
E = diag(V'*H*V);
end