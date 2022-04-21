function [ObjV_af,ObjV_cf,ObjV_ef,ObjV_sf,ObjV_tf]=find_ObjV(u,gen,ObjV_k1,ObjV_k2,...
         ObjV_a,ObjV_c,ObjV_e,ObjV_s,ObjV_t,ObjV_a2,ObjV_c2,ObjV_e2,ObjV_s2,ObjV_t2);      %求最优代的各项适应度
%gen=1时，父代为最初的父代；gen>1时，父代为重插后的父代
    [R1,T1]=find(u==ObjV_kgen);
    [R2,T2]=find(u==ObjV_k2gen);

    [R1,T1]=ismember(u,ObjV_kgen(:,gen));
    [R2,T2]=ismember(u,ObjV_k2);

    
if  R1~=0
    ObjV_af(gen,1)=ObjV_a(T1,1);
    ObjV_cf(gen,1)=ObjV_c(T1,1);
    ObjV_ef(gen,1)=ObjV_e(T1,1);
    ObjV_sf(gen,1)=ObjV_s(T1,1);
    ObjV_tf(gen,1)=ObjV_t(T1,1);
else
    ObjV_af(gen,1)=ObjV_a2(T2,1);
    ObjV_cf(gen,1)=ObjV_c2(T2,1);
    ObjV_ef(gen,1)=ObjV_e2(T2,1);
    ObjV_sf(gen,1)=ObjV_s2(T2,1);
    ObjV_tf(gen,1)=ObjV_t2(T2,1);
end

    
    
    