%�Ӻ������г��ܵ�����ʱ��
function [T_t,t]=func_ttotal(a,hs,v,NIND,NAVR);
t=zeros(NIND,NAVR);
for i=1:NIND
    for j=1:NAVR-1
        if  a(i,j)==0
            if  v(i,j)==0
                v(i,j)=1;       %����鵽a=0���������˶�ʱv=0����ʱ��Ϊ����v=1m/s
            end
            t(i,j)=abs((hs(i,j+1)-hs(i,j))/v(i,j));
        else
            t(i,j)=abs((sqrt(abs(v(i,j)^2+2*a(i,j)*(hs(i,j+1)-hs(i,j)))-v(i,j)))/a(i,j));
        end
    end
end
T_t=sum(t,2);                   %������ʱ��



