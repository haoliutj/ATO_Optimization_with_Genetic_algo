%�Ӻ�����������Ӧ�Ⱥ���
function ObjV_c=func_c(v,NIND,NAVR);
v_1=22.22;v_2=16.67;        %��·���٣���վ����
cc=zeros(NIND,NAVR-1);
for i=1:NIND
    for j=1:NAVR-2
        if  v(i,j)<v_1
            cc(i,j)=0;
        else
            cc(i,j)=v(i,j)-v_1;
        end
    end
    if  v(i,NAVR-1)<v_2
        cc(i,NAVR-1)=0;
    else
        cc(i,NAVR-1)=v(i,NAVR-1)-v_2;
    end
end
ObjV_c=sum(cc,2);

