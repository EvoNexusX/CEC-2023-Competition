function [peak,allpeak] = KEDE_v4(fn,run)
    if nargin == 0
        fn = 7;
        run= 5;
    end
    fprintf("正在进行Function %d run %d\n",fn,run);
    % 确定随机种子参数
    algRand = RandStream.create('mt19937ar','seed', run);
    RandStream.setGlobalStream(algRand);
    % 初始化问题参数
    pro = DMMOP(fn);
    % 初始化参数设置
    D = pro.D;
    min_popsize = 7 + floor(3*log(D));
    lambda = min_popsize;
    %  全局最优解集合
    bestmem_set = [];       
    bestval_set = [];
%     % 局部最优解集合
    subbestmem_set = [];    
    subbestval_set = [];  
    temp_best_pop  = [];

    first = 1;
    %% 判断问题是否终止
    while ~pro.Terminate()
         
        %% 基于密度生成种群，并评估适应值
        path = sprintf('./IDBPI_pop/init_pop_dim%d_run%d.mat',D,run);
        load(path);
        pop_I = pop;
        fits_I = pro.GetFits(pop_I);
        init_pop = [pop_I];
        fits = [fits_I];
        ss = get_ss(pro);
        if pro.D == 5
            if ~isempty(bestmem_set)
                old_pop = [];
                for i = 1:size(bestmem_set,1)
                    temp = bestmem_set(i,:)+ss*randn(algRand,lambda,pro.D);
                    old_pop = [old_pop;temp;bestmem_set(i,:)];
                end
                old_fit = pro.GetFits(old_pop);
                init_pop = [init_pop;old_pop];
                fits = [fits;old_fit];
            end
        end
        pop_F = FDBPI(pro.lower, pro.upper, init_pop, fits, 0.1*pro.freq, D,algRand);
        fits_F= pro.GetFits(pop_F); 
        init_pop = [init_pop;pop_F];
        fits = [fits; fits_F];
        if size(bestmem_set,1)~=1
            bestmem_set = bestmem_set(2:end,:);
        else
            bestmem_set = [];
        end
        %% 分簇
        [fits, sort_index] = sort(fits, 'descend');
        init_pop = init_pop(sort_index, :);
        species = NBC_lp(init_pop);
        num_species = length(species);
        species_arr = [species.len];
        imp_spec_count = length(species_arr(species_arr>=min_popsize));% 重要种群的数量
        
        %% 生成groups,并根据问题的复杂程度，初始化参数
        groups = struct();
        i=1;
        if first == 1
            if imp_spec_count>=9 % 说明问题比较复杂
                option = 1;
            else
                option = 0;
            end
        end
        if option == 0
            disp("----------------option 为  0 ------------------");
            for j = 1 : num_species

                groups(i).idx = i;
                groups(i).OPTS.first = 1;


                if species(j).len >= min_popsize
                    groups(i).OPTS.pop = init_pop(species(j).idx(1 : min_popsize), :);
                    groups(i).OPTS.val = fits(species(j).idx(1 : min_popsize));
                    groups(i).xmean = mean(groups(i).OPTS.pop)';
                    x = groups(i).OPTS.pop - groups(i).xmean';
                    groups(i).OPTS.sigma = sqrt((1/(min_popsize*D))*sum(x(:).^2));
                    groups(i).cc = std(groups(i).OPTS.val);
                    groups(i).bestval = fits(species(j).seed);
                    groups(i).bestmem = init_pop(species(j).seed, :);
                    groups(i).delta = 0;
                    groups(i).iters = 0;
                else
                    groups(i).OPTS.pop = init_pop(species(j).idx(1 : species(i).len), :);
                    groups(i).OPTS.val = fits(species(j).idx(1 : species(i).len));

                    if species(j).len == 1
                        groups(i).xmean = groups(i).OPTS.pop';
                        groups(i).OPTS.sigma = 0.5;
                    else
                        groups(i).xmean = mean(groups(i).OPTS.pop)'; %均值
                        x = groups(i).OPTS.pop - groups(i).xmean';
                        groups(i).OPTS.sigma = sqrt((1/((species(i).len)*D))*sum(x(:).^2)); % 方差
                    end

                    add_size = min_popsize - species(i).len;
                    sigma = groups(i).OPTS.sigma;

                    add_pop = groups(i).xmean' + sigma .* normrnd(0, 1, add_size, D);
                    add_fit = pro.GetFits(add_pop);

                    groups(i).OPTS.pop = [groups(i).OPTS.pop; add_pop];
                    groups(i).OPTS.val = [groups(i).OPTS.val; add_fit];
                    groups(i).cc = std(groups(i).OPTS.val);
                    groups(i).bestval = fits(species(j).seed);
                    groups(i).bestmem = init_pop(species(j).seed, :);
                    groups(i).delta = 0;
                    groups(i).iters = 0;
                end
                i = i + 1;
            end
            if ~isempty(bestmem_set)
                i  = i-1;
                for k = i+1:i+size(bestmem_set,1)
                    groups(k).idx = k;
                    groups(k).OPTS.first = 1;
                    groups(k).xmean = bestmem_set(k-i,:)';
                    groups(k).OPTS.pop = groups(k).xmean' + randn(algRand,lambda,pro.D);
                    x = groups(k).OPTS.pop - groups(k).xmean';
                    groups(k).OPTS.sigma = 0.01;
                    groups(k).bestmem = groups(k).xmean';
                    groups(k).bestval = 0;
                    groups(k).OPTS.val = groups(k).bestval;
                    groups(k).cc = 0.01;
                    groups(k).delta = 0;
                    groups(k).iters = 0;
                end
            end
        else
            for j = 1 : num_species
                if species(j).len < min_popsize
                    continue;
                end
                groups(i).idx = i;
                groups(i).OPTS.first = 1;
                groups(i).OPTS.pop = init_pop(species(j).idx(1 : min_popsize), :);
                groups(i).OPTS.val = fits(species(j).idx(1 : min_popsize));
                groups(i).xmean = mean(groups(i).OPTS.pop)';
                x = groups(i).OPTS.pop - groups(i).xmean';
                groups(i).OPTS.sigma = sqrt((1/(min_popsize*D))*sum(x(:).^2));
                groups(i).cc = std(groups(i).OPTS.val);
                groups(i).bestval = fits(species(j).seed);
                groups(i).bestmem = init_pop(species(j).seed, :);
                groups(i).delta = 0;
                groups(i).iters = 0;

                i = i + 1;
            end
          
        end

  
        bestmem_set = zeros(1,D);
        bestval_set = [-999];
      
        if pro.D>5
            option = 0;
        end

 
        %% 版本更新处，根据适应值排序group
        if option == 0
            [~,index] = sort([groups.bestval],'descend');
            groups = groups(index);
            num_groups = length(groups);
            lambda = min_popsize;
            temp_pop = zeros(num_groups*lambda,pro.D);
            % 建立总种群
            for j =1:num_groups
                temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
            end
        end
        num_groups = length(groups);
        %% 初始化每个簇的贡献
        itermax = ceil((0.25*pro.freq)/(num_groups*min_popsize));
%         itermax = 20;
        for i = 1:num_groups
            if option == 1
                [new_groups] = CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax,algRand);
                groups(i) = new_groups;
            else
                [new_groups] = KE_CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax, algRand,temp_pop,i,temp_best_pop);
                groups(i) = new_groups;
                temp_pop((i-1)*lambda+1:i*lambda,:) = groups(i).OPTS.pop;
            end  
        end
        val = [groups.bestval];  pop = cat(1, groups.bestmem);
        [bestval, ibest] = max(val);  bestmem = pop(ibest, :);
        bestval = round(bestval);
        if option == 0
            [~,index] = sort([groups.bestval],'descend');
            groups = groups(index);% 根据种群最优解进行排序
            for j =1:num_groups
                temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
            end
            itermax = 20;
        end
       tag  = 1;
        %% 循环迭代
        while ~pro.CheckChange(bestmem_set,bestval_set)
            fprintf("评估次数:%d\n",pro.evaluated);
            val = [groups.bestval]; pop = cat(1, groups.bestmem);
            [~, first_idx] = max(val);
    
            delta = [groups.delta];
            expected_gen = ceil((bestval - val)./(delta./itermax));
            expected_gen(first_idx) = Inf;
            if ~isempty(bestmem_set)
                gdis = pdist2(pop, bestmem_set);
                gdis = min(gdis, [], 2);
                temp_arr = [expected_gen', -gdis];
                [~, idx] = sortrows(temp_arr);
            else
                randnum = randperm(length(groups));
                temp_arr = [expected_gen', randnum'];
                [~, idx] = sortrows(temp_arr);
            end
            
            second_idx = idx(1);
            itermax = 20;
            % 演化最好的一个子群
            i = first_idx;
 
            rest = pro.freq - rem(pro.evaluated, pro.freq);
            if pro.change == 1
                continue;
            end
            if min_popsize*itermax <= rest
                if option == 1
                    [new_groups] = CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax,algRand);
                else
                    [new_groups] = KE_CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax, algRand,temp_pop,i,temp_best_pop);
                end
            else
                itermax = min(itermax, floor(rest/min_popsize));
                if itermax>0
                    if option == 1
                        [new_groups] = CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax,algRand);
                    else
                        [new_groups] = KE_CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax, algRand,temp_pop,i,temp_best_pop);

                    end
                end
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                useless_pop = rand(algRand,rest, D) .* (pro.upper - pro.lower) + pro.lower;
                useless_fit = pro.GetFits(useless_pop);
                continue;
            end    
            % 将参数拷贝到结构体中
            groups(i) = new_groups;
        
            if groups(i).bestval > bestval
                bestmem = groups(i).bestmem;
                bestval = groups(i).bestval;
            end
            if option == 0
                num_groups = length(groups);
                [~,index] = sort([groups.bestval],'descend');
                groups = groups(index);

                lambda = min_popsize;
               
                % 建立总种群
                for j =1:num_groups
                    temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
                end

            end
            % 演化潜力第二大的子群
            i = second_idx;
            itermax = 20;
            rest = pro.freq - rem(pro.evaluated, pro.freq);
            if pro.change == 1
                continue;
            end
            if min_popsize*itermax <= rest
                if option == 1
                    [new_groups] = CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax,algRand);
                else
                    [new_groups] = KE_CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax, algRand,temp_pop,i,temp_best_pop);
                end
            else
                itermax = min(itermax, floor(rest/min_popsize));
        
                if itermax>0
                    if option == 1
                        [new_groups] = CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax,algRand);
                    else
                        [new_groups] = KE_CMA_ES(pro, groups(i), pro.lower, pro.upper, itermax, algRand,temp_pop,i,temp_best_pop);
                    end
                end
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
                useless_fit = pro.GetFits(useless_pop);
                continue;
            end   
             % 将参数拷贝到结构体中
            groups(i) = new_groups;
            if option == 0
                num_groups = length(groups);
                [~,index] = sort([groups.bestval],'descend');
                groups = groups(index);
                % 建立总种群
                for j =1:num_groups
                    temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
                end
                first_idx = 1;
            end
            if ((groups(first_idx).cc <= 1e-1||groups(first_idx).OPTS.sigma<=0.1) && bestval-groups(first_idx).bestval>2)
                if pro.D == 5
                    groups(first_idx).xmean = -pro.D/2+0.5*randn(algRand,D,1);
                else
                    groups(first_idx).xmean = -pro.D/2+0.5*randn(algRand,D,1);
                end
                disp("该组不收敛");
                %                itermax = itermax + 80;
                tag = 1;
                groups(first_idx).OPTS.sigma = 1*(pro.D/5)^2;
                % %                groups(1).xmean = randn(algRand,pro.D,1);
                %                     groups(1).xmean = (mean(bestmem_set,1)/(size(bestmem_set,1)-1))';
                groups(first_idx).OPTS.pop = groups(1).xmean' + randn(algRand,lambda,pro.D);

                groups(first_idx).cc = 1.5;
                %                groups(2:end) = [];
                groups(first_idx).iters = 0;
                num_groups = num_groups - 1;
                %                 groups(1).cc = 0.5;
                %                 disp(groups(1).xmean);

                %                 num_groups = num_groups-1;

            end
               
            %尽可能减小废物种群的迭代次数



            % 收集收敛的种群
            if groups(first_idx).cc<=1e-7 && abs(groups(first_idx).bestval - bestval) <= 1e-5
                p = 1;
                %                 if abs(groups(first_idx).bestval - bestval) <= 1e-3
                if ~isempty(bestmem_set)
                    dis_arr = pdist2(groups(first_idx).bestmem, bestmem_set);
                    if min(dis_arr) < 1e-3
                        %                                                         groups(first_idx) = [];
                        p = 0;
                        %                             continue;
                    end
                end
                if p == 1
                    bestmem_set = [bestmem_set; groups(first_idx).bestmem];
                    bestval_set = [bestval_set; groups(first_idx).bestval];

                end
                if option == 0
                    temp_best_pop = [temp_best_pop;groups(1).OPTS.pop];
                end

                groups(first_idx) = [];
                if isempty(groups)
                    rest = pro.freq - rem(pro.evaluated, pro.freq);
                    useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
                    useless_fit = pro.GetFits(useless_pop);
                    continue;
                end
                continue
                %                 end
            end
            if ~isempty(bestmem_set)
                dis_arr = pdist2(groups(first_idx).bestmem, bestmem_set);
                if min(dis_arr) < 1e-3
                    groups(first_idx) = [];
                end
            end
            if isempty(groups)
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
                useless_fit = pro.GetFits(useless_pop);
                continue;
            end
        end
        disp(ss);
        peak = length(bestval_set(abs(bestval_set-max(bestval_set))<1e-5));
        fprintf("---Function%d Run%d Env%d/60一共找到%d个峰---\n",fn,run,pro.env,peak);
%         if ~isempty(bestmem_set)
%             bestmem_set = [];
%             bestval_set = [];
%         else
%             bestmem_set = zeros(1,D);
%             bestval_set = [0];
%         end
        first = 0; 
        temp_best_pop = [];
    end
    [peak, allpeak] = pro.GetPeak();
    PR = sum(peak, 2) ./ sum(allpeak, 2);
    fprintf("1e-3精度下:")
    disp(peak(1,:));
    fprintf("1e-4精度下:")
    disp(peak(2,:));
    fprintf("1e-5精度下:")
    disp(peak(3,:));
    disp(allpeak);
    disp(PR);
    filename = sprintf('./PEAKS/KEDE_v4_F%d_runs%d.mat', fn, run);
    save(filename,'peak','allpeak',"PR");
end

 