%子函数：求舒适度适应度函数
function ObjV_a=func_a(a);
aa=diff(a,1,2);  %g的各行中前后列项之差;考虑a1~=0;a10=0;舒适度的加速度差的数量为NAVR-1；
aa=abs(aa);
ObjV_a=sum(aa,2); %求每条染色体加速度变化之和(舒适度）;每行求和！


