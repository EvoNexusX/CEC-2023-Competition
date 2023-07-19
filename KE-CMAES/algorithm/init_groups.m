function [groups,split] = init_groups(pro,lambda,init_pop,species,groups,fits,sort_index,first,imp_spec_count,bestmem_set,algRand)
%     if first ==1
%         if imp_spec_count>=9 % 说明问题比较复杂
%             option = 1;
%         else
%             option = 0;
%         end
%     end
    split = 1;
    i = 1;
    D = pro.D;
    
    min_popsize = lambda;
    for k = 1:length(sort_index)
        j = sort_index(k);
        if i>9
            break;
        end
        if species(j).len >= min_popsize
            groups(i).idx = j;
            groups(i).OPTS.first = 1;
            groups(i).OPTS.pop = init_pop(species(j).idx(1 : min_popsize), :);
            groups(i).OPTS.val = pro.GetFits(init_pop(species(j).idx(1 : min_popsize), :));
            groups(i).xmean = mean(groups(i).OPTS.pop)';
            x = groups(i).OPTS.pop - groups(i).xmean';
            groups(i).OPTS.sigma = sqrt((1/(min_popsize*D))*sum(x(:).^2));
            groups(i).OPTS.count = 0;
            groups(i).cc = std(groups(i).OPTS.val);
            groups(i).bestmem = init_pop(species(j).seed, :);
            groups(i).bestval = pro.GetFits(groups(i).bestmem);
            groups(i).delta = 0;
            groups(i).iters = 0;
        else
            if i > 9 && split==1
                break;
            else
                split = 0;
                groups(i).idx = j;
                groups(i).OPTS.first = 1;
                groups(i).OPTS.pop = init_pop(species(j).idx(1 : species(j).len), :);
                groups(i).OPTS.val = fits(species(j).idx(1 : species(j).len));
                groups(i).OPTS.count = 0;
                if species(j).len == 1
                    groups(i).xmean = groups(i).OPTS.pop';
                    groups(i).OPTS.sigma = 0.5;
                else
                    groups(i).xmean = mean(groups(i).OPTS.pop)'; %均值
                    x = groups(i).OPTS.pop - groups(i).xmean';
                    groups(i).OPTS.sigma = sqrt((1/((species(j).len)*D))*sum(x(:).^2)); % 方差
                end
    
                add_size = min_popsize - species(j).len;
                sigma = groups(i).OPTS.sigma;
    
                add_pop = groups(i).xmean' + sigma .* normrnd(0, 1, add_size, D);
                add_fit = pro.GetFits(add_pop);
    
                groups(i).OPTS.pop = [groups(i).OPTS.pop; add_pop];
                groups(i).OPTS.val = [groups(i).OPTS.val; add_fit];
                groups(i).cc = std(groups(i).OPTS.val);
                groups(i).bestmem = init_pop(species(j).seed, :);                
                groups(i).bestval = pro.GetFits(groups(i).bestmem);
                groups(i).delta = 0;
                groups(i).iters = 0;
                if ~isempty(bestmem_set)
                    i  = i-1;
                    for k = i+1:i+size(bestmem_set,1)
                        groups(k).idx = k;
                        groups(k).OPTS.first = 1;
                        groups(k).xmean = bestmem_set(k-i,:)';
                        groups(k).OPTS.pop = groups(k).xmean' + randn(algRand,lambda,pro.D);
                        x = groups(k).OPTS.pop - groups(k).xmean';
                        groups(k).OPTS.sigma = 0.5;
                        groups(k).OPTS.count = 0;
                        groups(k).bestmem = groups(k).xmean';
                        groups(k).bestval = 0;
                        groups(k).OPTS.val = groups(k).bestval;
                        groups(k).cc = 0.01;
                        groups(k).delta = 0;
                        groups(k).iters = 0;
                    end
                end
            end
        end
        i = i + 1;
    end
    num_groups = length(groups);
    for i = 1:num_groups
        for j = 1:num_groups
            distance = 0;
            if i~=j
                distance = distance + pdist2(groups(i).xmean',groups(j).xmean');
            end
        groups(i).mean_distance = distance/(num_groups-1);
    end
    
end