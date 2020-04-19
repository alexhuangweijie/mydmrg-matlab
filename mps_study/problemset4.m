clear
%%%%% Problem set 4 solutions:
% (a) define tensor and normalize
d = 6;
H = zeros(d,d,d,d,d);
for ni = 1:d
    for nj = 1:d
        for nk = 1:d
            for nl = 1:d
                for nm = 1:d
                    H(ni,nj,nk,nl,nm) = sqrt(ni+2*nj+3*nk+4*nl+5*nm);
                end
            end
        end
    end
end
nmH = norm(H(:));
H1 = H/nmH;

% (b) multi-stage decomposition
chi = 6;
[utemp,stemp,vtemp] = svd(reshape(H1,[d^2,d^3]),0);
A = reshape(utemp(:,1:chi),[d,d,chi]);
Htemp = reshape(stemp(1:chi,1:chi)*(vtemp(:,1:chi)'),[chi,d,d,d]);
[utemp,stemp,vtemp] = svd(reshape(Htemp,[chi*d,d^2]),0);
B = reshape(utemp(:,1:chi),[chi,d,chi]);
C = reshape(stemp(1:chi,1:chi)*(vtemp(:,1:chi)'),[chi,d,d]);
% check accuracy
H2 = ncon({A,B,C},{[-1,-2,1],[1,-3,2],[2,-4,-5]});
TrErr = norm(H1(:)-H2(:))