

clear

m1 = 2;
m2 = 3;n2 = 3;
m3 = 4;
H = rand(m2,n2);
C = rand(m1*n2*m3,1);
A = rand(m1*n2*m3,1);
H = H'+H;
counter = 1;
Index = zeros(m1*m2*m3*n2,2);
Index2 = zeros(m1*m2*m3*n2,4);
for i1 = 1:m1
    for i2 = 1:m2
        for i3 = 1:m3
            for j2 = 1:n2
                index1 = (i1-1)*m2*m3 + (i2-1)*m3 + i3;
                index2 = (i1-1)*n2*m3 + (j2-1)*m3 + i3;
                Index(counter,:) = [index1,index2];
                Index2(counter,:) = [i1,i2,i3,j2];
                counter = counter +1;  
            end
        end
    end
end
Indicator = zeros(m1*m2*m3*n2,2);
Indicator2 = zeros(m1*m2*m3*n2,4);





parfor i = 1:m1*m2*m3*n2
%     ii1 = Index(i,1);
%     ii2 = Index(i,2);
%     ii3 = Index(i,3);
%     jj2 = Index(i,4);
%     index1 = ii1*ii2*ii3;
%     index2 = ii1*jj2*ii3;
%     index1 = (ii1-1)*m2*m3 + (ii2-1)*m3 + ii3;
%     index2 = (ii1-1)*n2*m3 + (jj2-1)*m3 + ii3;
%     Indicator(i,:) = [index1,index2];
%     Indicator2(i,:) = [ii1,ii2,ii3,jj2];
      iii1 = Index(i,1);
      iii2 = Index(i,2);
      jj1 = Index2(i,4)
      ii2 = Index2(i,2)
      A(iii1,iii2) = H(ii2,ji2)*C(iii2,1);
end
Indicator = unique(Indicator,'rows');
Indicator2 = unique(Indicator2,'rows');



