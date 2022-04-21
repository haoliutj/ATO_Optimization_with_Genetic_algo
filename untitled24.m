clear all
%�����Ŵ��㷨����
NIND=10;        %������Ŀ(Number of individuals)
MAXGEN=500;     %����Ŵ�����(Maximum number of generations)
NVAR=1;         %��������
PRECI=3;       %�����Ķ�����λ��(Precision of variables)
GGAP=0.9;       %����(Generation gap)
trace=zeros(2, MAXGEN);                        %Ѱ�Ž���ĳ�ʼֵ
FieldD=[3;-1;1;1;0;0;0];                    %����������(Build field descriptor)
Chrom=crtbp(NIND, PRECI);                      %��ʼ��Ⱥ
gen=0;                                         %��������
vari_mj=bs2rv(Chrom, FieldD);                 %�����ʼ��Ⱥ��ʮ����ת��
NN=length(vari_mj);
aa=[vari_mj(1);vari_mj(2);vari_mj(3);vari_mj(4);vari_mj(5);vari_mj(6);vari_mj(7);vari_mj(8);vari_mj(9);vari_mj(10)];
a=zeros(NN,1);
bb=[-1/9;-5/9;-1/3;-7/9;7/9;5/9;1/9;1/3;7/9;-1/3];
cc=[0.4;0.2;0.6;0;0.8;-0.25;-0.45;-0.65;0.8;0.6];
jj=ones(NN,1);
kk=ones(NN,1);
[jj(1),kk(1)]=find(bb==aa(1));
