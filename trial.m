clear
clc
format long
t=1;
noftk=40;%留下的向量的最大数目
niter=29;%链的长度
Bz=0;j=1;deltaT = 0.2;%磁场
measure1 = zeros(niter*2, 10);%测量结果储存
measure2 = zeros(niter*2, 20);
%measure3 = zeros(niter*2, 10);
%block从2个site开始 循环一次系统加一个site superblock一共2*Niter+2个site
%参数
tot_sites=2*niter+2;
half_sites=niter+1;
sites_in_system=niter;
Energymin=zeros(1000,1);
Energy = 0;
ExactEnergy = -2;
%Sz=zeros(tot_sites,1);
%数据记录
%初始化数据   stuct H包含有六个数据 包含算符和哈密顿量(I,Sz,Sp,SM,H)以及一个数（basis_size)
BLOCK.basis_size = 2;
BLOCK.I = eye(2);
BLOCK.C = [0,1;0,0];
BLOCK.Cdag = [0,0;1,0];
BLOCK.H = [0,0;0,0];
BLOCK.N = [0,0;0,1];
EMPTY.basis_size = 1;
EMPTY.I = 1;
EMPTY.C = 1;
EMPTY.Cdag = 1;
EMPTY.H = 1;
EMPTY.N = 1;
%参数初始化  把每个格点的初始参数输入到sys_block和environment_block中进行储存
%在长链的首和尾增加一个单位矩阵方便后面处理 所以对于链里面第i个格点 存在sys_block{i+1}里面
system_block = cell( tot_sites + 2 , 2 );environment_block = cell( tot_sites + 2 , 2 );
for i_1 = 1 : tot_sites+2
    %BLOCK.H=[0,0;0,-10*(mod(i_1,2)-0.5)];
    system_block{ i_1 , 1 } = BLOCK;environment_block{ i_1 , 1 } = BLOCK;
end
for i_2 = 1 : tot_sites+2
    system_block{ i_2 , 2 } = BLOCK;environment_block{ i_2 , 2 } = BLOCK;
end
system_block{ 1 , 1 } = EMPTY;environment_block{ 1 , 1 } = EMPTY;
system_block{ end , 1 } = EMPTY;environment_block{ end , 1 } = EMPTY;
%添加扰动部分
K.basis_size = 2;
K.I = eye(2);
K.C = [0,1;0,0];
K.Cdag = [0,0;1,0];
K.H = [0,0;0,-20];
K.N = [0,0;0,1];
system_block{1+half_sites,2} = K;
system_block{-2+half_sites,2} = K;
system_block{1+half_sites,1} = K;
system_block{-2+half_sites,1} = K;
groundstate = cell(tot_sites+2,1);
transform_L = cell(tot_sites+2,1);transform_L{2,1}=eye(2);transform_L{1,1}=1;
transform_R = cell(tot_sites+2,1);transform_R{2,1}=eye(2);transform_R{1,1}=1;
%有限系统循环
Block_L = BLOCK;Block_R = BLOCK;
for i = 1 : niter
    %第i轮循环的时候 此时的Free_site为编号为1+i , niter*2-i-1 的两个格点
    Free_site_L = system_block{ i + 2 , 2 };
    Free_site_R = system_block{ niter * 2 + 3 - i , 2 };
    %对角化哈密顿量
    LastEnergy = Energy;
    [Psi,Energy] = lanczos2(Block_L , Free_site_L , Free_site_R , Block_R ,j);
    %函数是将左右块和两个自由格点当作输入参数，返回的是系统哈密顿量和两个并入自由格点的块
    [ H_super , Block_L , Block_R ] = BlockBuilding_hubbard ( Block_L , Free_site_L , Free_site_R , Block_R , j  , Bz );
    %[ Psi , Energy ] = eigs ( H_super , 1 , 'SR' );
    EnergyPerBond = ( Energy - LastEnergy ) / 2;
    %此函数是进行基变换，将左右两个块的哈密顿量和算符进行变换  输入的是基态 函数里面含有求密度矩阵和约化密度矩阵的操作
    [ Block_L , Block_R ] = Rotate_hubbard( Block_L , Block_R , Psi  , noftk );
    system_block{ i+2 , 1 } = Block_L;
    environment_block{ i+2 , 1 } = Block_R;
end

%有限部分热身
Energy2 = Energy;
sites_in_environment = tot_sites - sites_in_system - 2;
counter = 0;Mark2=[];
for sweep = 1 : 1
    for iw = 1 : 3%ii指标反映对应的扫的时候的移动方向 奇数向右 偶数向左
        change = [ 1 , -1 , 1  ];
        while sites_in_system ~= -1 && sites_in_environment ~= -1
            Block_L = system_block{ sites_in_system + 1 , 1 };
            Block_R = environment_block{sites_in_environment + 1 , 1 };
            Free_site_L = system_block{ sites_in_system + 2 , 2 };
            Free_site_R = system_block{ sites_in_system + 3 , 2 };
            %对角化superblock哈密顿量 求基态
            LastEnergy2 = Energy2;
            [Psi,Energy] = lanczos2(Block_L , Free_site_L , Free_site_R , Block_R ,j );
            %构造哈密顿量
            [ H_super , Block_L , Block_R ] = BlockBuilding_hubbard ( Block_L , Free_site_L , Free_site_R , Block_R , j ,Bz );
            %[ Psi , Energy2 ] = eigs( H_super , 1 , 'SR' );
            %储存基态
            groundstate{sites_in_system+1} = Psi;
            convergence = ( Energy2 - LastEnergy2 ) / 2;
            EnergyPerBond2 = Energy2 / ( tot_sites - 1 );
            %测量
            [ measure_temp1 , measure_temp2 ] = Measure_hubbard ( Block_L , Free_site_L , Free_site_R , Block_R , Psi );
            %求密度矩阵 约化密度矩阵 求变换矩阵 返回变换矩阵
            transform_input1 = transform_L{sites_in_system+1};
            transform_input2 = transform_R{sites_in_environment+1};
            [ Block_L,Block_R,transform1,transform2] = Rotate_hubbard_transformation( Block_L,Block_R,Psi,noftk);
            %存储变换矩阵
            if mod(iw,2)==1
                transform_L{sites_in_system+2} = transform1;
            else
                transform_R{sites_in_environment+2} = transform2;
            end
            %在每个sweep中 选择ii等于2的时候进行的测量结果储存
            if  iw == 1|| iw == 3
                system_block{ sites_in_system + 2 , 1 } = Block_L;
            elseif iw == 2
                environment_block{ sites_in_environment + 2 , 1 } = Block_R;
                measure1(sites_in_system+2, sweep) = measure_temp1;
            end
            %记录每个块里面的格点数变化  ii奇数change为1 对应于左边system扩大 右边environment缩小
            sites_in_environment = sites_in_environment - change(iw);
            sites_in_system = sites_in_system + change(iw);
            Mark2=[Mark2; sites_in_system,sites_in_environment,sweep,iw];
            counter = counter + 1 ;
            if iw == 3 && sites_in_system == half_sites-1
                Mark2=[Mark2;0,0,0,0];
                break
            end
        end
        if iw==1
            sites_in_environment = 1;
            sites_in_system = tot_sites-3;
        elseif iw==2
            sites_in_system = 1;
            sites_in_environment = tot_sites-3;
        end
    end
end

%有限部分循环 使用上一轮的波函数加速
Energy2 = Energy;
sites_in_environment = tot_sites - sites_in_system - 2;
counter = 0;Mark2=[];
Psi_lanczos = rand(noftk^2*4,1);
for sweep = 1 : 10
    for ii = 1 : 3%ii指标反映对应的扫的时候的移动方向 奇数向右 偶数向左
        change = [ 1 , -1 , 1  ];
        while sites_in_system ~= 1 && sites_in_environment ~= 1
            Block_L = system_block{ sites_in_system + 1 , 1 };
            Block_R = environment_block{sites_in_environment + 1 , 1 };
            Free_site_L = system_block{ sites_in_system + 2 , 2 };
            Free_site_R = system_block{ sites_in_system + 3 , 2 };
            %对角化superblock哈密顿量 求基态
            LastEnergy2 = Energy2;
            [Psi,Energy] = lanczos_Psi(Block_L , Free_site_L , Free_site_R , Block_R ,j , Psi_lanczos);
            %构造哈密顿量
            [ H_super , Block_L , Block_R ] = BlockBuilding_hubbard ( Block_L , Free_site_L , Free_site_R , Block_R , j ,Bz );
            %[ Psi , Energy2 ] = eigs( H_super , 1 , 'SR' );
            %储存基态
            groundstate{sites_in_system+1} = Psi;
            convergence = ( Energy2 - LastEnergy2 ) / 2;
            EnergyPerBond2 = Energy2 / ( tot_sites - 1 );
            %测量
            [ measure_temp1 , measure_temp2 ] = Measure_hubbard ( Block_L , Free_site_L , Free_site_R , Block_R , Psi );
            %求密度矩阵 约化密度矩阵 求变换矩阵 返回变换矩阵
            transform_input1 = transform_L{sites_in_system+1};
            transform_input2 = transform_R{sites_in_environment+1};
            [Block_L,Block_R,transform1,transform2,Psi_L,Psi_R] = Rotate_hubbard_transformation_fast(Block_L,Block_R,Psi,noftk,transform_input1,transform_input2);
            if  ii == 1||ii==3
                Psi_lanczos = Psi_L;
            elseif ii == 2||ii==4
                Psi_lanczos = Psi_R;
            end
            %存储变换矩阵
            if mod(ii,2)==1
                transform_L{sites_in_system+2} = transform1;
            else
                transform_R{sites_in_environment+2} = transform2;
            end
            %在每个sweep中 选择ii等于2的时候进行的测量结果储存
            if  ii == 1|| ii == 3
                system_block{ sites_in_system + 2 , 1 } = Block_L;
            elseif ii == 2
                environment_block{ sites_in_environment + 2 , 1 } = Block_R;
                measure1(sites_in_system+2, sweep) = measure_temp1;
            end
            %记录每个块里面的格点数变化  ii奇数change为1 对应于左边system扩大 右边environment缩小
            sites_in_environment = sites_in_environment - change(ii);
            sites_in_system = sites_in_system + change(ii);
            Mark2=[Mark2; sites_in_system,sites_in_environment,sweep,ii];
            counter = counter + 1 ;
            if ii == 3 && sites_in_system == half_sites-1
                Mark2=[Mark2;0,0,0,0];
                break
            end
        end
        if ii==1
            sites_in_environment = 1;
            sites_in_system = tot_sites-3;
        elseif ii==2
            sites_in_system = 1;
            sites_in_environment = tot_sites-3;
        end
    end
end
%时间演化部分

Energy3 = Energy;
%选择从左边块里面没有格点为时间演化的开始
sites_in_system = 0  ;
sites_in_environment = tot_sites - sites_in_system - 2;
counter2 = 0;Mark1=[];
%读取初始的态 和自由格点的算符
Free_site_L = system_block{ 8 , 2 };
Free_site_R = system_block{ 8 , 2 };
Psi = groundstate{sites_in_system+1};
for sweep = 1 : 20
    %iii表示扫的时候移动的方向
    for iii = 1 : 4
        change = [ 1 , -1 , 1 , -1 , 1 ];
        while counter2<3000
            %在扫到缩小的块里面没有格点的时候就换方向
            if (sites_in_system==0&&iii==2)||(sites_in_environment==0&&iii==1)||(sites_in_environment==0&&iii==3)
                Mark1=[Mark1;0,0,0,0,3,iii];
                break
            end
            %结束一次sweep
            if iii == 4 && sites_in_system == 0
                Mark1=[Mark1;0,0,0,0,2,iii];
                break
            end
            Block_L = system_block{ sites_in_system + 1 , 1 };
            Block_R = environment_block{sites_in_environment + 1 , 1 };
            %计算时间演化算符 iii==1和2进行时间演化  iii=3的时候进行测量
            %iii=4不操作（为了使得每个sweep的起点相同）
            if iii==1||iii==2
                U = Uoperator_hubbard( Block_L , Free_site_L , Free_site_R , Block_R , j , deltaT );
            else
                U = 1;
            end
            Psi = U * Psi;
            %测量
            if iii == 3
                [ measure_temp1 , measure_temp2 ] = Measure_hubbard_time ( Block_L , Free_site_L , Free_site_R , Block_R , Psi );
                measure2(sites_in_system+2, sweep) = measure_temp1;
            end
            %读取变换矩阵
            transform_input1 = transform_L{sites_in_system+1};
            transform_input2 = transform_R{sites_in_environment+1};
            %输入演化之后的态 进行矩阵变换 变换矩阵
            %Psi_L Psi_R左边扩大和右扩大的情况的态
            [ Psi_L , Psi_R ,transform1 , transform2 ] = Rotate_state( Block_L , Block_R , Psi  , noftk ,transform_input1 ,transform_input2 );
            %将扩大的块的变换矩阵进行储存
            if iii==1||iii==2
                if mod(iii,2)==1
                    transform_L{sites_in_system+2} = transform1;
                else
                    transform_R{sites_in_environment+2} = transform2;
                end
            end
            %更新态Psi
            if  iii == 1||iii==3
                Psi = Psi_L;
            elseif iii == 2||iii==4
                Psi= Psi_R;
            end
            sites_in_environment = sites_in_environment - change(iii);
            sites_in_system = sites_in_system + change(iii);
            counter2 = counter2 + 1;
            Mark1=[Mark1; sites_in_system,sites_in_environment,size(Psi),1,iii];
        end
    end
end



fprintf('%.8f\t%.8f\t%.8f\n%.8f\t%.8f\t%.8f\n%.8f\n', Energy, EnergyPerBond,ExactEnergy-EnergyPerBond,Energy2,EnergyPerBond2,ExactEnergy-EnergyPerBond2 ,convergence);
figure(1)
plot(measure1(2:end,10));
hold on
plot(measure2(2:end,1:13))
% %figure(2)
% for index = 1:13
%     figure(index)
% plot(measure2(2:end,index));
% end
