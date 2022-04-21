
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
variable=bs2rv(Chrom, FieldD);                 %计算初始种群的十进制转换
NN=length(variable);
aa=[variable(1);variable(2);variable(3);variable(4);variable(5);variable(6);variable(7);variable(8);variable(9);variable(10)];
a=zeros(NN,1);
for ii=1:NN
    if aa(ii)==-0.7778
        a(ii)=0;
    else if aa(ii)==-0.5556
            a(ii)=0.2;
        else if aa(ii)==-0.1111
                a(ii)=0.4;
            else if aa(ii)==-0.3333
                    a(ii)=0.6;
                else if aa(ii)==0.7778
                a(ii)=0.8;
                    else if aa(ii)==0.5556
                a(ii)=-0.25;
                        end
                    end
                end
            end
        end
    end
end
a(0)=0;
i=1:10;
ObjV=abs(a(i+1)-a(i));        %计算目标函数值
while gen<MAXGEN
   FitnV=ranking(ObjV);                                  %分配适应度值(Assign fitness values)         
   SelCh=select('sus', Chrom, FitnV, GGAP);               %选择
   SelCh=recombin('xovsp', SelCh, 0.7);                   %重组
   SelCh=mut(SelCh);                                      %变异
   variable=bs2rv(SelCh, FieldD);                         %子代个体的十进制转换
%    ObjVSel=abs(variable-1000);                            %计算子代的目标函数值
   i=1:9;
   ObjVSel=abs(a(i+1)-a(i));
   [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel); %重插入子代的新种群
   variable=bs2rv(Chrom, FieldD);
   gen=gen+1;                                             %代计数器增加
   %输出最优解及其序号，并在目标函数图像中标出，Y为最优解,I为种群的序号
   [Y, I]=min(ObjV);  
   hold on;
   %figure(1);
   plot(variable(I), Y, 'bo');
   trace(1,gen)=min(ObjV);                               %遗传算法性能跟踪
   trace(2,gen)=sum(ObjV)/length(ObjV);
end
variable=bs2rv(Chrom, FieldD);                            %最优个体的十进制转换
 hold on, grid;
 figure(2);
 plot(variable,ObjV,'b*');
figure(3);
plot(trace(1,:)); grid
 hold on;
 plot(trace(2,:),'-.');grid
legend('解的变化','种群均值的变化')