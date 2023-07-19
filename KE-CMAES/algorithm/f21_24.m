max_run = 15;
% for dim = [5 10]
%    delete(gcp('nocreate'));
%    parpool('local',max_run);
%    spmd(max_run)
%        siglerun_IDBPI(dim, labindex);
%    end
% end

for func = 21 : 24
   delete(gcp('nocreate'));
   parpool('local',max_run);
   spmd(max_run)
       disp(func),disp(labindex);
       KEDE_v4(func, labindex);
   end
end
% for func = 21 : 24
%    delete(gcp('nocreate'));
%    parpool('local',max_run);
%    spmd(max_run)
%        KEDE_v4(func, labindex);
%    end
% end
%   