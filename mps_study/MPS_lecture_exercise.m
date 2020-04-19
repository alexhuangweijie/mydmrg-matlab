clear 
% bond dimensions
Da = 70; % alpha
Db = 71; % beta
Dc = 72; % gamma
Dd = 73; % delta
Dm = 74; % mu
A = rand(Dc,Dd); % tensor A(gamma,delta)
B = rand(Da,Dm,Dc); % tensor B(alpha,mu,gamma)
C = rand(Db,Dm,Dd); % tensor C(beta,mu,delta)
% tobj = tic2;
BC = contract(B,3,2,C,3,2); % BC(alpha,gamma,beta,delta)
ABC1 = contract(BC,4,[2 4],A,2,[1 2]); % ABC(alpha,beta)

ABC = ncon({A,B,C},{[1,2],[-1,3,1],[-2,3,2]})
sum(sum(ABC - ABC1))
