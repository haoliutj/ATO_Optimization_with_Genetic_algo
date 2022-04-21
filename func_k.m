%子函数：总的适应度函数
function ObjV_k=func_k(w_a,w_s,w_t,w_c,w_e,ObjV_a,ObjV_s,ObjV_t,ObjV_c,ObjV_e);
ObjV_k=w_a*ObjV_a+w_s*ObjV_s+w_t*ObjV_t+w_c*ObjV_c+w_e*ObjV_e;