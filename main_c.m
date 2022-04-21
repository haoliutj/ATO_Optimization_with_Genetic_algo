%��������
lbs=0;ubs=1000;   %�����Ա���x��Χ��0,1000��
lba=-1;uba=1;     %�����Ա���y��Χ��-1,1��
ps=10;pa=3;       %fieldD��ά��
m=200000;          %�г�����200��
w_c=0.4857;w_t=0.2728;w_s=0.1980;w_a=0.0303;w_e=0.0132; %����Ӧ�Ⱥ���Ȩ��ϵ��

%�����Ŵ��㷨����
NIND=100;         %������Ŀ(Number of individuals)
MAXGEN=20;        %����Ŵ�����(Maximum number of generations)
PRECI=13;         %�����Ķ�����λ��(Precision of variables)
NAVR=10;          %�ؼ�������
GGAP=0.7;         %����(Generation gap)
px=0.85;          %�������
pm=0.05;          %�������
Y=zeros(MAXGEN,1);I=zeros(MAXGEN,1); %�������Ž�Y�����Ŵ�I
trace=zeros(MAXGEN,2);                        %Ѱ�Ž���ĳ�ʼֵ
FieldD=[rep([ps pa],[1 NAVR]);rep([lbs lba],[1 NAVR]);rep([ubs uba],[1 NAVR]);...
    rep([1 1],[1,NAVR]);rep([0 0],[1,NAVR]);rep([1 0],[1 NAVR]);rep([1 0],[1 NAVR])];                    %����������(Build field descriptor)
Chrom=crtbp(NIND,PRECI*NAVR);                      %��ʼ��Ⱥ
%Chrom1=Chrom;
gen=0;                                         %��������

%����Ŀ�꺯���ӳ���1��������2����ѡ����̭����Ӵ���3�����ز����Ӵ�
variable1=bs2rv(Chrom, FieldD);                 %�����ʼ��Ⱥ��ʮ����ת��
hs=variable1(:,1:2:end); %λ��
ha=variable1(:,2:2:end); %��λ
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
ObjV_kgen=zeros(NIND,MAXGEN);%����ÿ��������ԭ�������ز��ĸ����ǣ��а�������������Ӧ�ȵľ���
ObjV_k2gen=[];                %zeros(GGAP*NIND,MAXGEN);  %����ɸѡ������а�������������Ӧ�ȵľ���

while gen<MAXGEN
   gen=gen+1;%������������
   ObjV_kgen(:,gen)=ObjV_k;
   FitnV=ranking(ObjV_k);                                  %������Ӧ��ֵ(Assign fitness values)         
   SelCh=select('rws', Chrom, FitnV, GGAP);               %ѡ��
   SelCh=recombin('xovdp', SelCh, px);                   %����
   SelCh=mut(SelCh,pm);                                      %����
   
   %ѡ����̭���Ӵ�Ŀ�꺯���ӳ���2����ѡ����̭����Ӵ�
   variable2=bs2rv(SelCh, FieldD);                         %�Ӵ������ʮ����ת��
   variable2=double(variable2);                             %���������Ͷ���Ϊdouble
   NIND2=length(variable2);
   hs2=variable2(:,1:2:end); %ѡ����̭���Ӵ�λ��
   ha2=variable2(:,2:2:end); %ѡ����̭���Ӵ���λ
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

   %�ز���Ӵ�Ŀ�꺯���ӳ���
   [Chrom ObjV_k]=reins(Chrom, SelCh, 1, 1, ObjV_k, ObjV_k2); %�ز���Ӵ�������Ⱥ
%    variable3=bs2rv(Chrom3, FieldD);
%    hs3=variable3(:,1:2:end); %�ز���Ӵ�λ��
%    ha3=variable3(:,2:2:end); %�ز���Ӵ���λ
%    ObjV_kgen=zeros(MAXGEN,1); 
%    ObjV_kgen(gen,1)=ObjV_k;
   %��ȡÿ�������Ž⼰����ţ�YΪ���Ž�,IΪ��������
   [u,w]=min(ObjV_k);                                    %����������u��w��u������¼A��ÿ�е���Сֵ��w��¼ÿ����Ⱥ����Ӧ����С�ĸ���λ��
   Chrom_optimize=Chrom(w,:);                             %��Ӧ�����ŵĸ���
   variable_f=bs2rv(Chrom_optimize, FieldD); %���Ÿ���ת��Ϊʮ���ƣ��õ����Ŵ����������λ���뵵λ
   hs_f=variable_f(:,1:2:end); %���Ŵ���λ��
   ha_f=variable_f(:,2:2:end); %���Ŵ��ĵ�λ
%    [ObjV_af,ObjV_cf,ObjV_ef,ObjV_sf,ObjV_tf]=find_ObjV(u,gen,ObjV_k1,ObjV_k2,...
%     ObjV_a,ObjV_c,ObjV_e,ObjV_s,ObjV_t,ObjV_a2,ObjV_c2,ObjV_e2,ObjV_s2,ObjV_t2);  %�����Ŵ��ĸ�����Ӧ��
  
   Y(gen,1)=u;
   I(gen,1)=w;
   trace(gen,1)=Y(gen,1);                            %����ÿ��������ֵ�Ĵ���
   trace(gen,2)=I(gen,1);                                      %����ÿ��������ֵ
end
% figure(1);
% plot(1:MAXGEN,trace(1:gen,1));
figure(3);
plot(trace(:,1));
