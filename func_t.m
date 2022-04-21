%子函数：求准时适应度函数
function ObjV_t=func_t(T_t,NIND);
T_target=90*ones(NIND,1);
ObjV_t=abs(T_t-T_target);

