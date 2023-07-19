max_run = 10;
delete(gcp('nocreate'));
parpool('local',max_run);
spmd(max_run)
    pro = DMMOP(17);
    D = pro.D;
    IDBPI(pro.lower, pro.upper, 0.4*pro.freq, D, labindex+10);
end
