function [psivec,dval] = My_eigs_2site(psivec,OPTS,linFunct,functArgs)
if norm(psivec) == 0
    psivec = rand(length(psivec),1);
end
psi = zeros(numel(psivec),OPTS.krydim+1);
A = zeros(OPTS.krydim,OPTS.krydim);
for k = 1:OPTS.maxit    
    psi(:,1) = psivec(:)/norm(psivec);
    for p = 2:OPTS.krydim+1
        psi(:,p) = linFunct(psi(:,p-1),functArgs{(1:length(functArgs))});
        for g = 1:1:p-1
            A(p-1,g) = dot(psi(:,p),psi(:,g));
            A(g,p-1) = conj(A(p-1,g));
        end
        for g = 1:1:p-1%与之前的向量保持正交化
            psi(:,p) = psi(:,p) - dot(psi(:,g),psi(:,p))*psi(:,g);
            psi(:,p) = psi(:,p)/max(norm(psi(:,p)),1e-16);
        end
    end
    
    [utemp,dtemp] = eig(0.5*(A+A'));
    xloc = find(diag(dtemp) == min(diag(dtemp)));
    psivec = psi(:,1:OPTS.krydim)*utemp(:,xloc(1));
end
psivec = psivec/norm(psivec);
dval = dtemp(xloc(1),xloc(1));
end