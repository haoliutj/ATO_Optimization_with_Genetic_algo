%基本参数
lbs=0;ubs=1000;   %函数自变量x范围【0,1000】
lba=-1;uba=1;     %函数自变量y范围【-1,1】
ps=10;pa=3;       %fieldD的维数
m=200000;          %列车重量200吨
w_c=0.4857;w_t=0.2728;w_s=0.1980;w_a=0.0303;w_e=0.0132; %总适应度函数权重系数

%定义遗传算法参数
NIND=100;         %个体数目(Number of individuals)
MAXGEN=20;        %最大遗传代数(Maximum number of generations)
PRECI=13;         %变量的二进制位数(Precision of variables)
NAVR=10;          %关键工况数
GGAP=0.7;         %代沟(Generation gap)
px=0.85;          %交叉概率
pm=0.05;          %变异概率
Y=zeros(MAXGEN,1);I=zeros(MAXGEN,1); %定义最优解Y与最优代I
trace=zeros(MAXGEN,2);                        %寻优结果的初始值
FieldD=[rep([ps pa],[1 NAVR]);rep([lbs lba],[1 NAVR]);rep([ubs uba],[1 NAVR]);...
    rep([1 1],[1,NAVR]);rep([0 0],[1,NAVR]);rep([1 0],[1 NAVR]);rep([1 0],[1 NAVR])];                    %区域描述器(Build field descriptor)
Chrom=crtbp(NIND,PRECI*NAVR);                      %初始种群
%Chrom1=Chrom;
gen=0;                                         %代计数器

%父代目标函数子程序；1代表父代；2代表选择淘汰后的子代；3代表重插后的子代
variable1=bs2rv(Chrom, FieldD);                 %计算初始种群的十进制转换
hs=variable1(:,1:2:end); %位置
ha=variable1(:,2:2:end); %档位
a=find_a(ha,NIND,NAVR);
ObjV_a=func_a(a);
ObjV_s=func_s(hs,NIND,NAVR);
v=func_v(a,hs,NIND,NAVR);
[T_t,t]=func_ttotal(a,hs,v,NIND,NAVR);
ObjV_t=func_t(T_t,NIND);
ObjV_c=func_c(v,NIND,NAVR);
ObjV_e=func_e(m,a,hs,NIND,NAVR);
ObjV_k=func_k(w_a,w_s,w_t,w_c,w_e,ObjV_a,ObjV_s,ObjV_t,ObjV_c,ObjV_e);
%ObjV_k1=ObjV_k;
ObjV_kgen=zeros(NIND,MAXGEN);%定义每个父代（原父代与重插后的父代们）中包含各个个体适应度的矩阵
ObjV_k2gen=[];                %zeros(GGAP*NIND,MAXGEN);  %定义筛选后各代中包含各个个体适应度的矩阵

while gen<MAXGEN
   gen=gen+1;%代计数器增加
   ObjV_kgen(:,gen)=ObjV_k;
   FitnV=ranking(ObjV_k);                                  %分配适应度值(Assign fitness values)         
   SelCh=select('rws', Chrom, FitnV, GGAP);               %选择
   SelCh=recombin('xovdp', SelCh, px);                   %重组
   SelCh=mut(SelCh,pm);                                      %变异
   
   %选择淘汰后子代目标函数子程序；2代表选择淘汰后的子代
   variable2=bs2rv(SelCh, FieldD);                         %子代个体的十进制转换
   variable2=double(variable2);                             %将数据类型定义为double
   NIND2=length(variable2);
   hs2=variable2(:,1:2:end); %选择淘汰后子代位置
   ha2=variable2(:,2:2:end); %选择淘汰后子代档位
   a2=find_a(ha2,NIND2,NAVR);
   ObjV_a2=func_a(a2);
   ObjV_s2=func_s(hs2,NIND2,NAVR);
   v2=func_v(a2,hs2,NIND2,NAVR);
   [T_t2,t2]=func_ttotal(a2,hs2,v2,NIND2,NAVR);
   ObjV_t2=func_t(T_t2,NIND2);
   ObjV_c2=func_c(v2,NIND2,NAVR);
   ObjV_e2=func_e(m,a2,hs2,NIND2,NAVR);
   ObjV_k2=func_k(w_a,w_s,w_t,w_c,w_e,ObjV_a2,ObjV_s2,ObjV_t2,ObjV_c2,ObjV_e2);
   ObjV_k2gen(:,gen)=ObjV_k2;

   %重插后子代目标函数子程序
   [Chrom ObjV_k]=reins(Chrom, SelCh, 1, 1, ObjV_k, ObjV_k2); %重插后子代的新种群
%    variable3=bs2rv(Chrom3, FieldD);
%    hs3=variable3(:,1:2:end); %重插后子代位置
%    ha3=variable3(:,2:2:end); %重插后子代档位
%    ObjV_kgen=zeros(MAXGEN,1); 
%    ObjV_kgen(gen,1)=ObjV_k;
   %获取每代的最优解及其序号，Y为最优解,I为个体的序号
   [u,w]=min(ObjV_k);                                    %返回行向量u和w，u向量记录A的每列的最小值，w记录每代种群中适应度最小的个体位置
   Chrom_optimize=Chrom(w,:);                             %适应度最优的个体
   variable_f=bs2rv(Chrom_optimize, FieldD); %最优个体转换为十进制，得到最优代各工况点的位置与档位
   hs_f=variable_f(:,1:2:end); %最优代的位置
   ha_f=variable_f(:,2:2:end); %最优代的档位
%    [ObjV_af,ObjV_cf,ObjV_ef,ObjV_sf,ObjV_tf]=find_ObjV(u,gen,ObjV_k1,ObjV_k2,...
%     ObjV_a,ObjV_c,ObjV_e,ObjV_s,ObjV_t,ObjV_a2,ObjV_c2,ObjV_e2,ObjV_s2,ObjV_t2);  %求最优代的各项适应度
  
   Y(gen,1)=u;
   I(gen,1)=w;
   trace(gen,1)=Y(gen,1);                            %记下每代的最优值的代数
   trace(gen,2)=I(gen,1);                                      %记下每代的最优值
end
% figure(1);
% plot(1:MAXGEN,trace(1:gen,1));
figure(3);
plot(trace(:,1));
