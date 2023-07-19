max_run = 10;
 

for func = 1 : 8
   delete(gcp('nocreate'));
   parpool('local',max_run);
   spmd(max_run)
       disp(func),disp(labindex);
       KEDE_v4(func, labindex);
   end
   delete(gcp('nocreate'));
end
% for func = 1 : 8
%    delete(gcp('nocreate'));
%    parpool('local',max_run);
%    spmd(max_run)
%        KEDE_v4(func, labindex+15);
%    end
% end
%   