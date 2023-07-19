max_run = 15;
% for dim = [5 10]
%    delete(gcp('nocreate'));
%    parpool('local',max_run);
%    spmd(max_run)
%        siglerun_IDBPI(dim, labindex);
%    end
% end

for func = 9 : 16
%    delete(gcp('nocreate'));
   parpool('local',max_run);
   spmd(max_run)
       disp(func),disp(labindex);
       KEDE_v4(func, labindex+15);
   end
end
% for func = 9 : 16
%    delete(gcp('nocreate'));
%    parpool('local',max_run);
%    spmd(max_run)
%        KEDE_v4(func, labindex+15);
%    end
% end
%   