function     [A] = mykron_parfor(H,m1,m2,m3,C)

% m1 = 2;
% m2 = 20;n2 = 20;
% m3 = 4;
% H = rand(m2,n2);
%C = rand(m1*n2*m3,1);
A = zeros(m1*m2*m3,1);
n2 = m2;
for i1 = 1:m1
    for i2 = 1:m2
        for  i3 = 1:m3
            parfor j2 = 1:m2
                index1 = (i1-1)*m2*m3 + (i2-1)*m3 + i3;
                index2 = (i1-1)*n2*m3 + (j2-1)*m3 + i3;
                A(index1,1) = A(index1,1)+H(i2,j2)*C(index2,1);
            end
        end
    end
end
% parfor i = 1:m1*m2*m3*n2-1
% i1 = floor(i/(m2*m3*n2))-1;
% z = [z ,i1];
% end


% A1 = kron(kron(eye(m1),H),eye(m3))*C;

 end
% delta = A-A1;
