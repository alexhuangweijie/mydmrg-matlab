function    [HL1,HL2,HR1,HR2] = nonhermitian(H1,H2)
HL1 = H1;
HL2 = H2;
HR1 = H1;
HR2 = H2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HL = H*H^T    get U
%HR = H^T*H    get V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HL1.H = H1.H*(H1.H).';
HL1.C = H1.C*(H1.C).';
HL1.Cdag = H1.Cdag*(H1.Cdag).';
HL1.N = H1.N*(H1.N).';
HL2.H = H2.H*(H2.H).';
HL2.C = H2.C*(H2.C).';
HL2.Cdag = H2.Cdag*(H2.Cdag).';
HL2.N = H2.N'*(H2.N).';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HR1.H = (H1.H).'*H1.H;
HR1.C = (H1.C).'*H1.C;
HR1.Cdag = (H1.Cdag).'*H1.Cdag;
HR1.N = (H1.N).'*H1.N;
HR2.H = (H2.H).'*H2.H;
HR2.C = (H2.C).'*H2.C;
HR2.Cdag = (H2.Cdag).'*H2.Cdag;
HR2.N = (H2.N).'*H2.N;
end