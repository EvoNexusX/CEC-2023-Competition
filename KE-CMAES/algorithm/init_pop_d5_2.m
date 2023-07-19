%%init_pop为IDBPI生成初始种群
max_run = 15;
delete(gcp('nocreate'));
parpool('local',max_run);
spmd(max_run)
    pro = DMMOP(1);% D为5
    D = pro.D;
    IDBPI(pro.lower, pro.upper, 0.4*pro.freq, D, labindex+15);
    fprintf("Run%d is init over!",labindex);
end 
