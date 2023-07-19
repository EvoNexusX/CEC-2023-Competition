func_nums=2;
func_runs=30;
BR = 0;
WR = 1;
for i = 22
    PRs = 0;
    for j = 1 : func_runs
        path='./PEAKS/KEDE_v4_F%d_runs%d.mat';
        file=sprintf(path, i , j );
        load(file);
%         disp(peak);
        if j == 15
            disp(PR);
        end
        PRs = PRs + PR;
        BR = max(BR,PR);
        WR = min(WR,PR);
    end
    PRs = PRs / func_runs;
    disp(PRs)
    disp(BR);
    disp(WR);
    
    pr_e1 = PRs(1,:);
    pr_e2 = PRs(2,:);
    pr_e3 = PRs(3,:);
    br_e1 = BR(1,:);
    br_e2 = BR(2,:);
    br_e3 = BR(3,:);
    wr_e1 = WR(1,:);
    wr_e2 = WR(2,:);
    wr_e3 = WR(3,:);

  
    range=strcat('A',num2str(i));
    writematrix(pr_e1,'./result/1e-3/pr2.xlsx',"Sheet",1,'Range',range);
    writematrix(pr_e2,'./result/1e-4/pr2.xlsx',"Sheet",1,'Range',range);
    writematrix(pr_e3,'./result/1e-5/pr2.xlsx',"Sheet",1,'Range',range);
    writematrix(br_e1,'./result/1e-3/br2.xlsx',"Sheet",1,'Range',range);
    writematrix(br_e2,'./result/1e-4/br2.xlsx',"Sheet",1,'Range',range);
    writematrix(br_e3,'./result/1e-5/br2.xlsx',"Sheet",1,'Range',range);
    writematrix(wr_e1,'./result/1e-3/wr2.xlsx',"Sheet",1,'Range',range);
    writematrix(wr_e2,'./result/1e-4/wr2.xlsx',"Sheet",1,'Range',range);
    writematrix(wr_e3,'./result/1e-5/wr2.xlsx',"Sheet",1,'Range',range);



end