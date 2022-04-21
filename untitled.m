
%定义遗传算法参数
NIND=10;        %个体数目(Number of individuals)
MAXGEN=30;     %最大遗传代数(Maximum number of generations)
NVAR=1;         %变量个数
PRECI=3;       %变量的二进制位数(Precision of variables)
GGAP=0.9;       %代沟(Generation gap)
trace=zeros(2, MAXGEN);                        %寻优结果的初始值
FieldD=[3;-1;1;1;0;0;0];                    %区域描述器(Build field descriptor)
Chrom=crtbp(NIND, PRECI);                      %初始种群
gen=0;                                         %代计数器
%计算父代的目标函数值
variable=bs2rv(Chrom, FieldD);                 %计算初始种群的十进制转换
variable=vpa(variable,4);                      %定义有效位数
NN=length(variable);
aa=variable;
a=zeros(NN,1);
bb=[-0.1111;-0.5556;-0.3333;-0.7778;0.7778;0.5556;0.1111;0.3333;0.7778;-0.3333];
cc=[0.4;0.2;0.6;0;0.8;-0.25;-0.45;-0.65;0.8;0.6];
[R,T]=ismember(aa,bb);
for ii=1:NN
    a(ii)=cc(T(ii),R(ii));
end
d=zeros(11,1);
for jj=2:11
    d(jj)=a(jj-1);
end
ObjV=diff(d);
ObjV=abs(ObjV);


while gen<MAXGEN
   FitnV=ranking(ObjV);                                  %分配适应度值(Assign fitness values)         
   SelCh=select('sus', Chrom, FitnV, GGAP);               %选择
   SelCh=recombin('xovsp', SelCh, 0.7);                   %重组
   SelCh=mut(SelCh);                                      %变异
   %计算子代的目标值
   variable=bs2rv(SelCh, FieldD);                         %子代个体的十进制转换
   variable=vpa(variable,4);                              %定义有效位数
   NN=length(variable);
   aa=variable;
   a=zeros(NN,1);
   bb=[-0.1111;-0.5556;-0.3333;-0.7778;0.7778;0.5556;0.1111;0.3333;0.7778;-0.3333];
   cc=[0.4;0.2;0.6;0;0.8;-0.25;-0.45;-0.65;0.8;0.6];
   [R,T]=ismember(aa,bb);
   for ii=1:NN
       a(ii)=cc(T(ii),R(ii));
   end
   
   d=zeros(10,1);
   for jj=2:10
       d(jj)=a(jj-1);
   end
   ObjVSel=diff(d);
   ObjVSel=abs(ObjV);
   ObjVSel=abs(diff(d));                                  %计算子代的目标函数值
   
   [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel); %重插入子代的新种群
   variable=bs2rv(Chrom, FieldD);
   gen=gen+1;                                             %代计数器增加
   %获取每代的最优解及其序号，Y为最优解,I为个体的序号
   [Y, I]=min(ObjV);                                    %返回行向量Y和I，Y向量记录A的每列的最小值，I向量记录每列最小值的行号
   trace(1,gen)=variable(I);                            %记下每代的最优值的代数
   trace(2,gen)=Y;                                      %记下每代的最优值
end
figure(1);
plot(trace(1,:),trace(2,:),'bo');                       %画出每代的最优点
grid on;
xlabel('遗传代数')
ylabel('最优解的变化')
title('进化过程')
plot(variable,ObjV,'b*');   %画出最后一代的种群
hold off

variable=bs2rv(Chrom, FieldD);                          %最优个体的十进制转换
 hold on, grid;
 figure(2);
 plot(variable,ObjV,'b*');
figure(3);
plot(trace(1,:)); grid
 hold on;
 plot(trace(2,:),'-.');grid
legend('解的变化','种群均值的变化')

[Y,I]=min(ObjV);
   trace(1,gen)=X(I);                            %记下每代的最优值
   trace(2,gen)=Y;                               %记下每代的最优值
end
plot(trace(1,:),trace(2,:),'bo');                            %画出每代的最优点
grid on;
plot(X,ObjV,'b*');   %画出最后一代的种群
hold off
%% 画进化图
figure(2);
plot(1:MAXGEN,trace(2,:));
grid on
xlabel('遗传代数')
ylabel('解的变化')
title('进化过程')
bestY=trace(2,end);
bestX=trace(1,end);
fprintf(['最优解:\nX=',num2str(bestX),'\nY=',num2str(bestY),'\n'])