function [pop] = IDBPI(lb, ub, init_popsize, dim, runs)
    %UNTITLED 此处显示有关此函数的摘要
    %   此处显示详细说明
    algRand = RandStream.create('mt19937ar','seed', runs); % 将run设置为随机种子
    RandStream.setGlobalStream(algRand);
    pop = zeros(init_popsize, dim);
    
    Cd = (pi^(dim/2))/(gamma(dim/2 + 1));
    R = (((1/(init_popsize))*prod(ub - lb))/Cd)^(1/dim);
    
    h = R;
    k = 10;
    
    % 阈值
    rho_min = exp(-h^2 / (2*h^2));
    
    pop(1, :) = lb + (ub - lb) .* rand(algRand,1, dim);
    ppl = lb + (ub - lb) .* rand(algRand,k*(init_popsize-1) ,dim);      % 计划探测位置
    
    for i = 2 : init_popsize
        temp_location = ppl(k*(i-2)+1:k*(i-1), :);
        dist = pdist2(temp_location, pop(1:i-1, :));
        rho = exp(-dist.^2 ./ (2 * h^2));
        rho(rho < rho_min) = 0;
        rho = (1/(i-1)) * sum(rho, 2);
        [~, selected_idx] = min(rho);
        pop(i, :) = temp_location(selected_idx, :);
        
        if mod(i, 5000) == 0
            pro = i/init_popsize;
            fprintf('%2d.%2d| progress: %.2f  \n', dim, runs, pro);
        end
    end
    init_pop_path=sprintf("./IDBPI_pop2/init_pop_dim%d_run%d",dim,runs);
    save(init_pop_path,"pop");
end