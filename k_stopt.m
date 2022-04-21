
clear all
%�����Ŵ��㷨����
NIND=100;        %������Ŀ(Number of individuals)
MAXGEN=500;     %����Ŵ�����(Maximum number of generations)
NVAR=1;         %��������
PRECI=10;       %�����Ķ�����λ��(Precision of variables)
GGAP=0.9;       %����(Generation gap)
trace=zeros(2, MAXGEN);                        %Ѱ�Ž���ĳ�ʼֵ
FieldD=[10;0;1024;1;0;1;1];                    %����������(Build field descriptor)
Chrom=crtbp(NIND, PRECI);                      %��ʼ��Ⱥ
gen=0;                                         %��������
variable=bs2rv(Chrom, FieldD);                 %�����ʼ��Ⱥ��ʮ����ת��
ObjV=abs(variable-1000);        %����Ŀ�꺯��ֵ
while gen<MAXGEN
   FitnV=ranking(ObjV);                                  %������Ӧ��ֵ(Assign fitness values)         
   SelCh=select('sus', Chrom, FitnV, GGAP);               %ѡ��
   SelCh=recombin('xovsp', SelCh, 0.7);                   %����
   SelCh=mut(SelCh);                                      %����
   variable=bs2rv(SelCh, FieldD);                         %�Ӵ������ʮ����ת��
   ObjVSel=abs(variable-1000);             %�����Ӵ���Ŀ�꺯��ֵ
   [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel); %�ز����Ӵ�������Ⱥ
   variable=bs2rv(Chrom, FieldD);
   gen=gen+1;                                             %������������
   %������Ž⼰����ţ�����Ŀ�꺯��ͼ���б����YΪ���Ž�,IΪ��Ⱥ�����
   [Y, I]=min(ObjV);  
   hold on;
   %figure(1);
   plot(variable(I), Y, 'bo');
   trace(1,gen)=min(ObjV);                               %�Ŵ��㷨���ܸ���
   trace(2,gen)=sum(ObjV)/length(ObjV);
end
variable=bs2rv(Chrom, FieldD);                            %���Ÿ����ʮ����ת��
 hold on, grid;
 figure(2);
 plot(variable,ObjV,'b*');
figure(3);
plot(trace(1,:)); grid
 hold on;
 plot(trace(2,:),'-.');grid
legend('��ı仯','��Ⱥ��ֵ�ı仯')