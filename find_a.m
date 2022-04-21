% 子函数：寻找加速度
function a=find_a(ha,NIND,NAVR);
ha=vpa(ha,4);                      %定义有效位数
a=zeros(NIND,NAVR);
bb=[-0.1111;-0.5556;-0.3333;-0.7778;0.7778;0.5556;0.1111;0.3333];
cc=[0.4;0.2;0.6;0;0.8;-0.25;-0.45;-0.65];
[R,T]=ismember(ha,bb);
for ii=1:NAVR
    for jj=1:NIND
        a(jj,ii)=cc(T(jj,ii),R(jj,ii));
    end
end


