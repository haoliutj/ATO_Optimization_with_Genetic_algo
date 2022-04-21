%子函数：列车总的运行时间
function [T_t,t]=func_ttotal(a,hs,v,NIND,NAVR);
t=zeros(NIND,NAVR);
for i=1:NIND
    for j=1:NAVR-1
        if  a(i,j)==0
            if  v(i,j)==0
                v(i,j)=1;       %当检查到a=0，即匀速运动时v=0，此时人为给定v=1m/s
            end
            t(i,j)=abs((hs(i,j+1)-hs(i,j))/v(i,j));
        else
            t(i,j)=abs((sqrt(abs(v(i,j)^2+2*a(i,j)*(hs(i,j+1)-hs(i,j)))-v(i,j)))/a(i,j));
        end
    end
end
T_t=sum(t,2);                   %运行总时间



