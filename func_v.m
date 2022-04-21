%子函数：求速度
function v=func_v(a,hs,NIND,NAVR);
v=zeros(NIND,NAVR);
for jj=1:NIND
    for ii=2:NAVR
        v(jj,ii)=sqrt(abs(2*a(jj,ii-1)*(hs(jj,ii)-hs(jj,ii-1))+v(jj,ii-1)^2));
    end
end

