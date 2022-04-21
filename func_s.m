%子函数：求停车精度适应度函数
function ObjV_s=func_s(hs,NIND,NAVR);
s_target=1000*ones(NIND,NAVR);
ObjV_s=abs(hs(:,NAVR)-s_target(:,NAVR));

