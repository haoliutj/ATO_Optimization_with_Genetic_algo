%�Ӻ������ܺ���Ӧ�Ⱥ���
function ObjV_e=func_e(m,a,hs,NIND,NAVR);
ee=zeros(NIND,NAVR-1);
for i=1:NIND
    for j=2:NAVR
        ee(i,j-1)=m*a(i,j-1)*abs((hs(i,j)-hs(j-1)));
    end
end
ObjV_e=sum(ee,2);     %�Ը������
