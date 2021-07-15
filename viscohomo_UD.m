function viscohomo_UD(nelx,nely,volfrac,penal,rmin,ft,omega)
%2021/6/8 16:21
%regu 阶段变化，与beta相关,为了体分比约束控制
%ISO 去除分数乘子
%% INITIALIZES
f_11=10;%给tan11的权系数
f_22=0;%给tan22的权系数
ff=1;%给tan11/tan22的权系数
c_f=25;%刚度约束

xP=1; % 1=中心圆/2=矩形/3=给定图像/4=随机矩阵‘
etah=volfrac;

plane=1;%平面应力=1/平面应变=2
ISO=1; %实部各向同性=1/虚部和实部都是=2
epsilon=0.001; %等式约束的平方小于某个小量
eps_b=0.003;%beta 更新的约束参数
C_beta=1.2;
ar=0.5;%即初始为单元数目的1/10；
regu0=nely*ar;%为了使得体积分数约束放在目标中和目标同量级
b_loop=50;
lu=0.3;%low=x-lu.up=x+lu;

M_name=mfilename;%获取当前M文件的名称
root_directory = 'E:\Visco_homo\UD';
hms = clock;  % hour, minute and second
pathName = [root_directory '\' date '-' num2str(hms(4)) '-' num2str(hms(5)) '-' num2str(ceil(hms(6)))...
    '-' M_name  '-' num2str(f_11) '_' num2str(f_22) '_' num2str(ff) '_' num2str(ISO)  '\'];
dircommand=['mkdir ' pathName];
system(dircommand);%创建了存储的文件
copyfile([M_name, '.m'],pathName) %保存当前M文件

span=3;
phi=90;

%保存优化迭代输出记录
diary([pathName, 'log.txt'] );
diary on
%不变% MATERIAL PROPERTIES/材料参数
syms s t w
Et1=70;
v1=0.22;
Et2=1+2.5*exp(-t);
v2=0.35;

Es1=s*laplace(Et1,t,s);
Es2=s*laplace(Et2,t,s);
Ew1=subs(Es1,s,1i*w);
Ew2=subs(Es2,s,1i*w);
lx=50;
ly=50;
dx=lx/nelx;
dy=ly/nely;

E1=eval(subs(Ew1,w,omega));
E2=eval(subs(Ew2,w,omega));

%采用拉梅常数的变换公式得到lambda和mu.
Lambda1=[E1*v1/(1+v1)/(1-2*v1) E2*v2/(1+v2)/(1-2*v2)];
Mu=[E1/2/(1+v1) E2/2/(1+v2)];
if plane==1 %平面应力
    Lambda=2*Mu.*Lambda1./(Lambda1+2.*Mu);%改成平面应力的情况
elseif plane==2 %平面应变
    Lambda=Lambda1;
end
nel=nelx*nely;

% D1=E1/(1-v1^2)*[1 v1 0;v1 1 0; 0 0 (1-v1)/2];
% D2=E2/(1-v2^2)*[1 v2 0;v2 1 0; 0 0 (1-v2)/2];
fprintf('f_11; f_22;ff;regu;lu \n'); %需要修改各向同性约束等
disp([f_11,f_22,ff,regu0,lu]);
%% 输出计算的参数/包括对称性参数
fprintf('volf;C2_R;C3_R;C11;C22;C2_I;C3_I\n'); %需要修改各向同性约束等
%C2 为C11=C22/加上C_orth=0得到正方对称的/加上C3得到各向同性结构
%C_orth=0对于二维结构可以使用一个强制对称轴实现
%正方对称可以使用两个对称轴实现
%各向同性可以由正方对称加入C3约束

%%%%%%%%强制X或者Y左右|上下对称==》各向异性结果
% axisy=1;%1=x;2=y;3=xy;0=无对称约束/ 其中x,y 为矩阵类型的坐标系
axisy=2;%1=x;2=y;3=xy;0=无对称约束/ 其中x,y 为矩阵类型的坐标系

% %%%只加入刚度约束和体积约束的结构
% vec=[1,0,0,1,0]; %带有体积约束/1刚度约束的
% vec=[1,0,0,1,1]; %带有体积约束/2刚度约束的
% vec=[0,0,0,1,0]; %无体积约束/1刚度约束的
% vec=[0,0,0,1,1]; %无体5积约束/2刚度约束的
%%%使用等式约束得到square symmetry结构
vec=[1,1,0,1,0]; %带有体积约束/刚度约束的%因为C2=0所以只需要一个刚度约束
% vec=[0,1,0,1,0];  %无体积约束/刚度约束的
%vec=[0,1,0,0,0]; %无体积约束/无刚度约束的

%%%使用等式约束得到各向同性的结构
% vec=[1,1,1,1,0]; %带有体积约束/刚度约束的%因为C2=0所以只需要一个刚度约束
% vec=[0,1,1,1,0];  %无体积约束/刚度约束的
%vec=[0,1,1,0,0]; %无体积约束/无刚度约束的%这里恐怕会得到全黑的结构

%%%%%%%%强制XY对称==》正方对称的结构
% axisy=3;%1=x;2=y;3=xy;/ 其中x,y 为矩阵类型的坐标系

% %%%只加入刚度约束和体积约束的结构square symmetry结构
% vec=[1,0,0,1,0]; %带有体积约束/1刚度约束的
% vec=[1,0,0,1,1]; %带有体积约束/2刚度约束的
% vec=[0,0,0,1,0]; %无体积约束/1刚度约束的
% vec=[0,0,0,1,1]; %无体积约束/2刚度约束的

%%%使用等式约束得到各向同性的结构/加入C3
% vec=[1,0,1,1,0]; %带有体积约束/刚度约束的%因为C2=0所以只需要一个刚度约束
%vec=[0,0,1,1,0];  %无体积约束/刚度约束的
%vec=[0,0,1,0,0]; %无体积约束/无刚度约束的%这里恐怕会得到全黑的结构

if ISO==2%虚部也要各向同性的时候
    vec=[vec 1  1];
else
    vec=[vec  0  0];
end

disp(vec); %为了查看使用的约束
fprintf('nelx,nely,volfrac,penal,rmin,ft,omega,axisy,c_f,plane,ISO,eps_b \n'); %需要修改各向同性约束等
disp([nelx,nely,volfrac,penal,rmin,ft,omega,axisy,c_f,plane,ISO,eps_b]);

% % 将结果的CH保存为txt文档
%将迭代的结果和重要的信息分别写入txt文档
file_CH=[pathName, 'CH.txt']; %写入变成A(:),即成为一行 而不是3*3
file_ite=[pathName, 'iteration.txt'];
file_sens=[pathName, 'sensitivity.txt'];
%% FINITE DISCRETIZATION/有限元离散
% Node numbers and element degrees of freedom for full (not periodic) mesh
[keLambda,keMu,feLambda,feMu]=elementMatVec(dx/2,dy/2,phi);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);% 节点的编号 矩阵形式的
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);%自由度的编号从节点开始（2n-2,2n）
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
%自由度编号 单元形式的 仍是逆时针 从左下角开始
%% IMPOSE PERIODIC BOUNDARY CONDITIONS/周期边界
% Use original edofMat to index into list with the periodic dofs
nn=(nelx+1)*(nely+1);%Total number of nodes
nnP=nelx*nely;%Total number of unique nodes
nnPArray=reshape(1:nnP,nely,nelx);
% Extend with a mirror of the top border
nnPArray(end+1,:)=nnPArray(1,:);
% Extend with a mirror of the left border
nnPArray(:,end+1)=nnPArray(:,1);
% Make a vector into which we can index using edofMa;
dofVector=zeros(2*nn,1);
dofVector(1:2:end)=2*nnPArray(:)-1;
dofVector(2:2:end)=2*nnPArray(:);
edofMat=dofVector(edofMat);
ndof=2*nnP;%Number of dofs
% Indexing vectors
iK=kron(edofMat,ones(8,1))';% 这里相当于复制，就是行索引，每个重复8次。
jK=kron(edofMat,ones(1,8))';
%% PREPARE FILTER/过滤矩阵
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION/迭代初始化
% The displacement vectors corresponding to the unit strain cases
chi0=zeros(nel,8,3);
% The element displacements for the three unit strains
chi0_e=zeros(8,3);
% 这里的ratio不重要 是因为1：2这个是写在矩阵里面而不是前面的系数，关于加载是单位应力就不需要再进行比例相除
%比如说diag[2,1,1]，那么之后要除以系数矩阵[4,2,2;2,1,1;2,1,1]，这样得到的才是Ch，即材料参数。
ke=keMu+keLambda;%Here the exact ratio does not matter because it is reflected in the load vector
fe=feMu+feLambda;%Here the exact ratio does not matter because it is reflected in the load vector
chi0_e([3 5:end],:)=ke([3 5:end],[3 5:end])\fe([3 5:end],:);
% epsilon0_11=(1 0 0)
chi0(:,:,1)=kron(chi0_e(:,1)',ones(nel,1));
% epsilon0_22=(0 1 0)
chi0(:,:,2)=kron(chi0_e(:,2)',ones(nel,1));
% epsilon0_12=(0 0 1)
chi0(:,:,3)=kron(chi0_e(:,3)',ones(nel,1));
CH=zeros(3);
DCH=cell(3,3);%对于密度的灵敏度
cellVolume=lx*ly;
%% 初始化密度矩阵x以及过滤后的物理量xPhys
switch (xP)
    case 1
        %%% 一般的圆
        x=repmat(volfrac,nely,nelx);
        ref=min(nelx,nely)*1/4;%初始化弱相
        for i=1:nelx
            for j=1:nely
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2)<ref
                    x(j,i)=1;
                end
            end
        end
%         a1=sum(sum(x))/volfrac/nelx/nely;
%         x=x./a1;
    case 2
        %%%Init_ele
        x=Init_ele(nelx,volfrac);
    case 3
        %         %%%采用给定的图像作为初始构型
        %         I=imread('chiral2.png');
        %         II=imresize(I,[nely,nelx]);
        %         x=double(rgb2gray(II))./255;
        %
        %%% 读取黑白图像作为初始构型：n_holes/cells
        I=imread('cells6.png');
        II=imresize(I,[nely,nelx]);
        rr=double(max(II(:)));
        x=double(II)./rr;
    case 4
        %%% 随机矩阵作为初始构型
        x=rand(nely,nelx);
end
%有密度/灵敏度/Heaviside过滤
beta=1;
if ft==1||ft==2
    xPhys=x;
elseif ft==3
    xTilde=x;
    [xh,~]=Heaviside(etah,beta,xTilde(:));
    xPhys=reshape(xh,nely,nelx);
end
%% GEOMERTRIC SYMMETRY/几何对称
% 将整个单胞化成四份%按照一般的矩阵的排列，左上角的右手坐标系
L1=(1:nely/2)';%左上边
L2=(nely/2+1:nely)';%左下边
U1=(1:nelx/2)';%上左边
U2=(nelx/2+1:nelx)';%上右边

x1=xPhys(L1,U1);%左上部分
x2=xPhys(L1,U2);%右上一部分
x3=xPhys(L2,U1);%左下
x4=xPhys(L2,U2);%右下
%对称轴设置axis=?;%1=x;2=y;3=xy; 其中x,y 为矩阵类型的坐标系x为|，y为→，为对称轴
switch (axisy)
    case 1 %关于x对称，左右相等
        xPhys=[x1 flip(x1,2); x3 flip(x3,2)];
    case 2 %关于y对称，上下对称
        xPhys=[x1 x2; flip(x1) flip(x2)];
    case 3 % 只有1/4独立
        xPhys=[x1 flip(x1,2); flip(x1) flip(flip(x1,2))];
    case 0
        xPhys=[x1 x2; x3 x4];%没有几何对称要求
end

loopbeta=0;
change=1;
loop=0;

figure('visible','off')
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal;drawnow;
xlim([0.5,nelx+0.5]);xticks([]);yticks([]);
print(gcf,'-dpng',[pathName,'design configuration' num2str(0) '.png']);
close(gcf)
%% 优化迭代
regu=regu0;
while (change >= 0.0005)&&(loop<=1000)
    loop=loop+1;
    loopbeta=loopbeta+1;
    %%  ASSEMBLE STIFFNESS MATRIX/组装刚度矩阵
    % Material properties in the different elements
    %此处采用线性插值/但是如果加入惩罚不知道还能否得到类似的结论？
    lambda=Lambda(1)*xPhys.^penal+Lambda(2)*(1-xPhys.^penal);
    mu=Mu(1)*xPhys.^penal+Mu(2)*(1-xPhys.^penal);
    % The corresponding stiffness matrix entries
    sK=keLambda(:)*lambda(:).'+keMu(:)*mu(:).';
    K=sparse(iK(:),jK(:),sK(:),ndof,ndof);
    %% LOAD VECTORS AND SOLUTION/加载以及有限元解
    % Assembly three load cases corresponding to the three strain cases
    sF=feLambda(:)*lambda(:).'+feMu(:)*mu(:).';
    iF=repmat(edofMat',3,1);
    jF=[ones(8,nel);2*ones(8,nel);3*ones(8,nel)];
    F=sparse(iF(:),jF(:),sF(:),ndof,3);
    % Solve (remember to constrain one node)
    chi(3:ndof,:)=K(3:ndof,3:ndof)\F(3:ndof,:);
    %% HOMOHENIZATION/均匀化
    for i=1:3
        for j=1:3
            sumLambda=((chi0(:,:,i)-chi(edofMat+(i-1)*ndof))*keLambda).*...
                (chi0(:,:,j)-chi(edofMat+(j-1)*ndof));
            sumMu=((chi0(:,:,i)-chi(edofMat+(i-1)*ndof))*keMu).*...
                (chi0(:,:,j)-chi(edofMat+(j-1)*ndof));
            sumLambda=reshape(sum(sumLambda,2),nely,nelx);
            sumMu=reshape(sum(sumMu,2),nely,nelx);
            % Homogenized elasticity tensor
            CH(i,j)=1/cellVolume*sum(sum(lambda.*sumLambda+mu.*sumMu));
            DCH{i,j}=1/cellVolume*(penal*xPhys.^(penal-1).*(Lambda(1)-Lambda(2)).*sumLambda...
                +penal*xPhys.^(penal-1).*(Mu(1)-Mu(2)).*sumMu);
        end
    end
    disp('---Homogenized elasticity tensor---');
    disp(CH);
    dlmwrite(file_CH,CH','-append','delimiter',' '); %此处变成行向量了
    %% 目标函数以及灵敏度分析
    CH_r=real(CH);
    CH_i=imag(CH);
    
    tan11=CH_i(1,1)/CH_r(1,1);
    tan22=CH_i(2,2)/CH_r(2,2);
    
    DchI_11=imag(DCH{1,1});
    DchR_11=real(DCH{1,1});
    DchI_12=imag(DCH{1,2});
    DchR_12=real(DCH{1,2});
    DchI_22=imag(DCH{2,2});
    DchR_22=real(DCH{2,2});
    DchI_33=imag(DCH{3,3});
    DchR_33=real(DCH{3,3});
    %%%实际上是不可能约束虚部的，因为各向同性只针对stiffness，unidirectional damper才有意义
    %%%预备是加入对称性约束 手动只优化1/4
    %%%基础的使用C11作为目标函数
    % 优化目标为拉压模量损耗之差
    CHr_1=CH_r(1,1);
    CHi_1=CH_i(1,1);
    CHr_2=CH_r(2,2);
    CHi_2=CH_i(2,2);
    %更改目标为max{tan11/tan22+tan11};
    tan_O=ff*tan11/tan22+f_11*tan11-f_22*tan22;
    dtan11=(DchI_11.*CHr_1-DchR_11.*CHi_1)./(CHr_1)^2;
    dtan22=(DchI_22.*CHr_2-DchR_22.*CHi_2)./(CHr_2)^2;
    Dtan_O=ff*(1/tan22^2*(dtan11.*tan22-dtan22.*tan11))+f_11.*dtan11-f_22.*dtan22;
    
    C2=CH(1,1)-CH(2,2);
    C3=CH(1,1)+CH(2,2)-2*(CH(1,2)+2*CH(3,3));
    C2_R=real(C2);
    C3_R=real(C3);
    C2_I=imag(C2);
    C3_I=imag(C3);
    
    DC2_R=DchR_11-DchR_22;
    DC3_R=DchR_11+DchR_22-2*(DchR_12+2*DchR_33);
    DC2_I=DchI_11-DchI_22;
    DC3_I=DchI_11+DchI_22-2*(DchI_12+2*DchI_33);
    
    %直接采用实部的结果,Ci^2<=0计算,d{f^2}=2*f*df
    %% 使用axisy强制对称约束的话/如果为0的话就不进行限制
    %等式约束的施加/f(x)=0中的fx本身就是个平方式
    dv=ones(nely,nelx)./nelx./nely; %这里是df(x)中括号的部分
    vol_diff=sum(xPhys(:))/nelx/nely-volfrac;%这里是体积差
    %%%vec(1)=0的话就不约束体积分数
    fval_v=regu*vol_diff^2;
    dfdx_v=regu*2*vol_diff*dv;
    tan_D=tan_O-vec(1)*fval_v;%所以目标是最小化-tan+(a*f(x))^2;
    Dtan=Dtan_O-vec(1)*dfdx_v;%灵敏度为-Dtan+2*(a*f(x))*(a*df(x))
    
    switch axisy
        case 1 %关于x对称，左右相等
            Dtan=sym_X(Dtan,L1,L2,U1);
            DchR_11=sym_X(DchR_11,L1,L2,U1);
            DchR_22=sym_X(DchR_22,L1,L2,U1);
            DC3_R=sym_X(DC3_R,L1,L2,U1);
            DC2_R=sym_X(DC2_R,L1,L2,U1);
            DC3_I=sym_X(DC3_I,L1,L2,U1);
            DC2_I=sym_X(DC2_I,L1,L2,U1);
        case 2 %关于y对称，上下对称
            Dtan=sym_Y(Dtan,L1,U1,U2);
            DchR_11=sym_Y(DchR_11,L1,U1,U2);
            DchR_22=sym_Y(DchR_22,L1,U1,U2);
            DC3_R=sym_Y(DC3_R,L1,U1,U2);
            DC2_R=sym_Y(DC2_R,L1,U1,U2);
            DC3_I=sym_Y(DC3_I,L1,U1,U2);
            DC2_I=sym_Y(DC2_I,L1,U1,U2);
        case 3 % 只有1/4独立
            Dtan=sym_XY(Dtan,L1,U1);
            DchR_11=sym_XY(DchR_11,L1,U1);
            DchR_22=sym_XY(DchR_22,L1,U1);
            DC3_R=sym_XY(DC3_R,L1,U1);
            DC2_R=sym_XY(DC2_R,L1,U1);
            DC3_I=sym_XY(DC3_I,L1,U1);
            DC2_I=sym_XY(DC2_I,L1,U1);
    end
    %% FILTERING/MODIFICATION OF SENSITIVITIES 灵敏度的链式求解
    if ft == 1%敏度过滤
        Dtan(:) = H*(x(:).*Dtan(:))./Hs./max(1e-3,x(:));
    elseif ft == 2%密度过滤也是需要链式法则的啊
        Dtan(:) = H*(Dtan(:)./Hs);
        DchR_11(:)=H*(DchR_11(:)./Hs);
        DchR_22(:)=H*(DchR_22(:)./Hs);
        DC3_R(:)=H*(DC3_R(:)./Hs);
        DC2_R(:)=H*(DC2_R(:)./Hs);
        DC3_I(:)=H*(DC3_I(:)./Hs);
        DC2_I(:)=H*(DC2_I(:)./Hs);
    elseif ft==3%heaviside投影，应该是过滤（模糊化）之后的映射（清晰化）
        [~,dxh]=Heaviside(etah,beta,xTilde(:));%xTilde就是密度过滤后的结果
        dx=reshape(dxh,nely,nelx);
        Dtan(:)=H*(Dtan(:).*dx(:)./Hs);
        DchR_11(:)=H*(DchR_11(:).*dx(:)./Hs);
        DchR_22(:)=H*(DchR_22(:).*dx(:)./Hs);
        DC3_R(:)=H*(DC3_R(:).*dx(:)./Hs);
        DC2_R(:)=H*(DC2_R(:).*dx(:)./Hs);
        DC3_I(:)=H*(DC3_I(:).*dx(:)./Hs);
        DC2_I(:)=H*(DC2_I(:).*dx(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES/更新变量
    f0val=-tan_D;
    fval1=C2_R^2-epsilon;  %CH_r(1,1)=CH_r(2,2);
    fval2=C3_R^2-epsilon;  %C3_R=0;
    fval3=-CH_r(1,1)+c_f;% D11>=15;
    fval4=-CH_r(2,2)+c_f;% D22>=15;
    fval5=C2_I^2-epsilon;  %CH_i(1,1)=CH_i(2,2);
    fval6=C3_I^2-epsilon;  %C3_I=0;
    
    df0dx=-Dtan(:);%max tan\delta =minimize (-tan\delta)
    dfdx1=2*C2_R.*DC2_R(:)';
    dfdx2=2*C3_R.*DC3_R(:)';
    dfdx3=-DchR_11(:)';
    dfdx4=-DchR_22(:)';
    dfdx5=2*C2_I.*DC2_I(:)';
    dfdx6=2*C3_I.*DC3_I(:)';
    
    f_all=[fval1; fval2;fval3;fval4;fval5;fval6];
    iter=[mean(xPhys(:));C2;C3;CH_r(1,1);CH_r(2,2);CH_r(3,3)]';
    dlmwrite(file_ite,iter,'-append','delimiter',' ');
    df_all=[dfdx1;dfdx2;dfdx3;dfdx4;dfdx5;dfdx6];
    DD=[-Dtan_O(:)';df_all];
    senswrite=[max(DD,[],2)  min(DD,[],2) mean(DD,2) ]'; %包含目标函数在内的所有灵敏度
    dlmwrite(file_sens,senswrite,'-append','delimiter',' ');
    
    T=diag(vec(2:end));%因为体积约束放在了目标里面 所以只有后面的部分
    T(all(T==0,2),:) = [];
    fval=T*f_all;
    dfdx=T*df_all;
    
    xval=x(:);
    xmin=zeros(length(xval),1);
    xmax=ones(length(xval),1);
    %     low=xmin;
    %     upp=xmax;
    low=max(xmin,xval-lu);
    upp=min(xmax,xval+lu);
    epsimin = 0.0000001;
    n=length(xval);
    m=size(fval,1);
    eeem    = ones(m,1);
    zerom   = zeros(m,1);
    cc=1000*eeem;
    d= eeem;
    a0= 1;
    a= zerom;
    raa0= 0.01;
    raa= 0.01*eeem;
    if loop==1
        xold1=xval;
        xold2=xval;
    end
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
        gcmmasub(m,n,loop,epsimin,xval,xmin,xmax,low,upp, ...
        raa0,raa,f0val,df0dx,fval,dfdx,a0,a,cc,d); %#ok<ASGLU>
    
    %      [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
    %     mmasub(m,n,loop,xval,xmin,xmax,xold1,xold2, ...
    %     f0val,df0dx,fval,dfdx,low,upp,a0,a,cc,d);%#ok<ASGLU>
    xold1=xold2;
    xold2=xval;
    xnew=reshape(xmma,nely,nelx);
    %初始保证了对称性的x,在后续的优化中进行对称性的设置
    x11=xnew(L1,U1);%左上部分
    x22=xnew(L1,U2);%右上一部分
    x33=xnew(L2,U1);%左下
    x44=xnew(L2,U2);%右下
    
    switch (axisy)
        case 1 %关于x对称，左右相等
            xnew_sym=[x11 flip(x11,2); x33 flip(x33,2)];
        case 2 %关于y对称，上下对称
            xnew_sym=[x11 x22; flip(x11) flip(x22)];
        case 3 % 只有1/4独立
            xnew_sym=[x11 flip(x11,2); flip(x11) flip(flip(x11,2))];
        case 0
            xnew_sym=[x11 x22; x33 x44];%没有几何对称要求
    end
    
    if ft == 1
        xPhys = xnew_sym;
    elseif ft == 2
        xPhys(:) = (H*xnew_sym(:))./Hs;
    elseif ft==3
        xTilde(:)=(H*xnew_sym(:))./Hs;
        [xh,~]=Heaviside(etah,beta,xTilde(:));
        xPhys=reshape(xh,nely,nelx);
    end
    change = max(abs(xnew_sym(:)-xval(:)));
    x = xnew_sym;
    %% 强制对称性
    x1=xPhys(L1,U1);%左上部分
    x2=xPhys(L1,U2);%右上一部分
    x3=xPhys(L2,U1);%左下
    x4=xPhys(L2,U2);%右下
    %对称轴设置axis=?;%1=x;2=y;3=xy; 其中x,y 为矩阵类型的坐标系x为|，y为→，为对称轴
    switch (axisy)
        case 1 %关于x对称，左右相等
            xPhys=[x1 flip(x1,2); x3 flip(x3,2)];
        case 2 %关于y对称，上下对称
            xPhys=[x1 x2; flip(x1) flip(x2)];
        case 3 % 只有1/4独立
            xPhys=[x1 flip(x1,2); flip(x1) flip(flip(x1,2))];
        case 0
            xPhys=[x1 x2; x3 x4];%没有几何对称要求
    end
    %% PRINT RESULTS
    fprintf(' It.:%5i   tan11.:%11.4f   tan22.:%11.4f  C2_R.:%11.4f C3_R.:%11.4f  C2_I.:%11.4f C3_I.:%11.4f  Vol.:%7.3f  ch.:%7.3f\n',...
        loop,tan11,tan22,C2_R,C3_R,C2_I,C3_I,mean(xPhys(:)),change);
    %% PLOT DENSITIES
    figure('visible','off')
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; drawnow;
    xlim([0.5,nelx+0.5]);xticks([]);yticks([]);
    if loop==1
        cplot1=tan11;
        vplot=mean(xPhys(:));
        clow1=0;
        cup1=cplot1+20;
    else
        cplot1=[cplot1,tan11]; %#ok<*AGROW>
        clow1=min(cplot1(:));
        cup1=max(cplot1(:));
        vv=mean(xPhys(:));
        vplot=[vplot,vv];
    end
    
    if mod(loop,span)==0
        print(gcf,'-dpng',[pathName,'design configuration' num2str(loop/span) '.png']);
    end
    clf;
    
    figure('visible','off')
    yyaxis left
    plot(1:loop,cplot1');
    xlabel('Iteration');
    ylabel('Objective tan\delta_1_1');
    ylim([clow1 cup1]);
    yyaxis right
    plot(1:loop,vplot');
    ylabel('Volfraction');
    ylim([0 1]);
    legend('Compliance','Volfraction');
    print(gcf,'-dpng',[pathName,'Iteration.png']);
    close(gcf)
    
    if mod(loop,span)==0
        save([pathName,'xPhy' num2str(loop/span) '.mat'],'xPhys');
    end
    clf;
    if ft==3&&beta<100&&(loopbeta>=b_loop||change<=eps_b)
        beta=C_beta*beta;
        loopbeta=0;
        change=1;
        fprintf('Parameter beta/regu increased to %g/%g.\n',beta,regu);
    end
    if beta>10&&beta<20
        regu=5*regu0;
    end
    if beta>20&&beta<30
        regu=15*regu0;
    end
    if beta>30
        regu=30*regu0;
    end
    
end
diary off
return
end

