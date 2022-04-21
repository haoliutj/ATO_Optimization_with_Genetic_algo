
%�����Ŵ��㷨����
NIND=10;        %������Ŀ(Number of individuals)
MAXGEN=30;     %����Ŵ�����(Maximum number of generations)
NVAR=1;         %��������
PRECI=3;       %�����Ķ�����λ��(Precision of variables)
GGAP=0.9;       %����(Generation gap)
trace=zeros(2, MAXGEN);                        %Ѱ�Ž���ĳ�ʼֵ
FieldD=[3;-1;1;1;0;0;0];                    %����������(Build field descriptor)
Chrom=crtbp(NIND, PRECI);                      %��ʼ��Ⱥ
gen=0;                                         %��������
%���㸸����Ŀ�꺯��ֵ
variable=bs2rv(Chrom, FieldD);                 %�����ʼ��Ⱥ��ʮ����ת��
variable=vpa(variable,4);                      %������Чλ��
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
   FitnV=ranking(ObjV);                                  %������Ӧ��ֵ(Assign fitness values)         
   SelCh=select('sus', Chrom, FitnV, GGAP);               %ѡ��
   SelCh=recombin('xovsp', SelCh, 0.7);                   %����
   SelCh=mut(SelCh);                                      %����
   %�����Ӵ���Ŀ��ֵ
   variable=bs2rv(SelCh, FieldD);                         %�Ӵ������ʮ����ת��
   variable=vpa(variable,4);                              %������Чλ��
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
   ObjVSel=abs(diff(d));                                  %�����Ӵ���Ŀ�꺯��ֵ
   
   [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel); %�ز����Ӵ�������Ⱥ
   variable=bs2rv(Chrom, FieldD);
   gen=gen+1;                                             %������������
   %��ȡÿ�������Ž⼰����ţ�YΪ���Ž�,IΪ��������
   [Y, I]=min(ObjV);                                    %����������Y��I��Y������¼A��ÿ�е���Сֵ��I������¼ÿ����Сֵ���к�
   trace(1,gen)=variable(I);                            %����ÿ��������ֵ�Ĵ���
   trace(2,gen)=Y;                                      %����ÿ��������ֵ
end
figure(1);
plot(trace(1,:),trace(2,:),'bo');                       %����ÿ�������ŵ�
grid on;
xlabel('�Ŵ�����')
ylabel('���Ž�ı仯')
title('��������')
plot(variable,ObjV,'b*');   %�������һ������Ⱥ
hold off

variable=bs2rv(Chrom, FieldD);                          %���Ÿ����ʮ����ת��
 hold on, grid;
 figure(2);
 plot(variable,ObjV,'b*');
figure(3);
plot(trace(1,:)); grid
 hold on;
 plot(trace(2,:),'-.');grid
legend('��ı仯','��Ⱥ��ֵ�ı仯')

[Y,I]=min(ObjV);
   trace(1,gen)=X(I);                            %����ÿ��������ֵ
   trace(2,gen)=Y;                               %����ÿ��������ֵ
end
plot(trace(1,:),trace(2,:),'bo');                            %����ÿ�������ŵ�
grid on;
plot(X,ObjV,'b*');   %�������һ������Ⱥ
hold off
%% ������ͼ
figure(2);
plot(1:MAXGEN,trace(2,:));
grid on
xlabel('�Ŵ�����')
ylabel('��ı仯')
title('��������')
bestY=trace(2,end);
bestX=trace(1,end);
fprintf(['���Ž�:\nX=',num2str(bestX),'\nY=',num2str(bestY),'\n'])