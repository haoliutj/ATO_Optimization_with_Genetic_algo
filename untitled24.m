clear all
%定义遗传算法参数
NIND=10;        %个体数目(Number of individuals)
MAXGEN=500;     %最大遗传代数(Maximum number of generations)
NVAR=1;         %变量个数
PRECI=3;       %变量的二进制位数(Precision of variables)
GGAP=0.9;       %代沟(Generation gap)
trace=zeros(2, MAXGEN);                        %寻优结果的初始值
FieldD=[3;-1;1;1;0;0;0];                    %区域描述器(Build field descriptor)
Chrom=crtbp(NIND, PRECI);                      %初始种群
gen=0;                                         %代计数器
vari_mj=bs2rv(Chrom, FieldD);                 %计算初始种群的十进制转换
NN=length(vari_mj);
aa=[vari_mj(1);vari_mj(2);vari_mj(3);vari_mj(4);vari_mj(5);vari_mj(6);vari_mj(7);vari_mj(8);vari_mj(9);vari_mj(10)];
a=zeros(NN,1);
bb=[-1/9;-5/9;-1/3;-7/9;7/9;5/9;1/9;1/3;7/9;-1/3];
cc=[0.4;0.2;0.6;0;0.8;-0.25;-0.45;-0.65;0.8;0.6];
jj=ones(NN,1);
kk=ones(NN,1);
[jj(1),kk(1)]=find(bb==aa(1));
