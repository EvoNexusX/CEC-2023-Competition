function [group,track] = KE_CMA_ES(pro, group, lb, ub, itermax, algRand,temp_pop,idx,temp_best_pop,imp_spec_count,track)
    
    xmean   = group.xmean;
    bestmem = group.bestmem;
    bestval = group.bestval;
    OPTS    = group.OPTS;
%     idx     = group.idx;
    old_bestval = group.bestval;
    lambda = 7 + floor(3*log(pro.D));
%     if idx == 1
%         disp(size(temp_pop,1));
%         res_pop = temp_pop(lambda+1:size(temp_pop,1),:);
%     else
%         if idx == size(temp_pop,1)
%         res_pop = temp_pop(1:(idx-1)*lambda,:)+temp_pop(idx*lambda+1:size(temp_pop,1),:);
% 
%     end
    temp_pop((idx-1)*lambda+1:end,:)=[];
    temp_pop = [temp_pop;temp_best_pop];
    % 初始化系数
    dim = length(xmean);
    sigma = OPTS.sigma;
    
    if OPTS.first == 1
        lambda = 7 + floor(3*log(dim));
%         lambda = length(group);
        mu = floor(lambda/2);
        % Strategy parameter setting: Selection
        weights = log(mu+1/2)-log(1:mu)';       % muXone recombination weights
        mu = floor(mu);                         % number of parents/points for recombination
        weights = weights/sum(weights);         % normalize recombination weights array
        mueff=sum(weights)^2/sum(weights.^2);   % variance-effective size of mu
        
        % Strategy parameter setting: Adaptation
        cc = (4+mueff/dim) / (dim+4 + 2*mueff/dim);     % time constant for cumulation for C
        cs = (mueff+2)/(dim+mueff+5);                   % t-const for cumulation for sigma control
        c1 = 2 / ((dim+1.3)^2+mueff);                   % learning rate for rank-one update of C
        cmu = 2 * (mueff-2+1/mueff) / ((dim+2)^2+2*mueff/2);    % and for rank-mu update
        damps = 1 + 2*max(0, sqrt((mueff-1)/(dim+1))-1) + cs;   % damping for sigma
        
        % Initialize dynamic (internal) strategy parameters and constants
        pc = zeros(dim,1); ps = zeros(dim,1);           % evolution paths for C and sigma
        B = eye(dim);                                   % B defines the coordinate system
        D = eye(dim);                                   % diagonal matrix D defines the scaling
        C = B*D*(B*D)';                                 % covariance matrix
        chiN=dim^0.5*(1-1/(4*dim)+1/(21*dim^2));        % expectation of ||N(0,I)|| == norm(randn(N,1))
        countval = 0;
        iters = 0;
    else
       
        lambda = OPTS.lambda;
        weights = OPTS.weights;
        mu = OPTS.mu;
        mueff = OPTS.mueff;
        cc = OPTS.cc;
        cs = OPTS.cs;
        c1 = OPTS.c1;
        cmu = OPTS.cmu;
        damps = OPTS.damps;
        pc = OPTS.pc;
        ps = OPTS.ps;
        B = OPTS.B;
        D = OPTS.D;
        C = OPTS.C;
        chiN = OPTS.chiN;
        countval = OPTS.countval;
        iters = group.iters;
    end
    % -------------------- Generation Loop --------------------------------
    stopiters = iters + itermax; 
    % Fes = 0;
    while iters < stopiters
        % Generate and evaluate lambda offspring
        
%         if OPTS.first == 1
%             arx = OPTS.pop';
%             arfitness = OPTS.val';
%             for k = 1 : lambda
%                 arz(:, k) = pinv(D) * pinv(B) * ((arx(:, k) - xmean)/sigma);
%             end
%         else
            % 生成10*lambda个


            temp_x =zeros(dim, 10);
%             temp_arx = zeros(dim, lambda)';
            for j = 1 : lambda
                temp_arz = randn(algRand, dim, 10);                          % standard normally distributed vector
                for k = 1 : 10
                    temp_x(:,k) = xmean + sigma*(B*D*temp_arz(:,k));
                end
                rho = kernel(temp_x',temp_pop,dim,imp_spec_count);
                group.OPTS.count =group.OPTS.count + length(rho(rho>0));

                [~, index] = min(rho);
                arx(:, j) = temp_x(:, index);
                arz(:, j) = temp_arz(:,index);
            end
            if group.OPTS.count >= 15*lambda
                fprintf("----The group is useless!----!\r");
                group = [];
                break;
            end
%                 [~,index] = sort(rho);

            for k = 1: lambda
%                 arx(:,k) = xmean + sigma * (B*D * arz(:,k));
                
                
                temp_ub = arx(:, k) > ub(1);
                temp_lb = arx(:, k) < lb(1);
                if any(temp_ub) || any(temp_lb)
                    arx(temp_ub, k) =  ub(1) ;
                    arx(temp_lb, k) =  lb(1);
                    arz(:, k) = pinv(D) * pinv(B) * ((arx(:, k) - xmean)/sigma);
                end
                countval = countval + 1;
            end
            
            arfitness = pro.GetFits(arx');
            % Fes = Fes + lambda;
            iters = iters + 1;
%         end
        if size(arfitness,1)~=size(arx',1)
            break;
        end

        % Sort by fitness and compute weighted mean into xmean
        [arfitness, arindex] = sort(arfitness,'descend');         % maxmization
        xmean = arx(:,arindex(1:mu))*weights;
        track = [track;xmean'];
        zmean = arz(:,arindex(1:mu))*weights;
      
        % Cumulation: Update evolution paths
        ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*countval/lambda))/chiN < 1.4+2/(dim+1);
        pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B*D*zmean);
        
       

        % Adapt covariance matrix C
        C = (1-c1-cmu) * C ...
            + c1 * (pc*pc' ... % plus rank one update
            + (1-hsig) * cc*(2-cc) * C) ... % minor correction
            + cmu ... % plus rank mu update
            * (B*D*arz(:,arindex(1:mu))) ...
            * diag(weights) * (B*D*arz(:,arindex(1:mu)))';
        
        % Adapt step-size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        
        % Update B and D from C
        C = triu(C) + triu(C,1)';       % enforce symmetry
        [B,D] = eig(C);                 % eigen decomposition, B==normalized eigenvectors
        D = diag(sqrt(diag(abs(D))));   % D contains standard deviations now
        % Break, if fitness satisfies stop condition
        
        if OPTS.first == 1
            OPTS.first = 0;
        end
        
        if arfitness(1) > bestval
            bestmem = arx(:,arindex(1))';
            bestval = arfitness(1);
        end
        
        if std(arfitness) < 1e-7
            break;
        end
    end
    if ~isempty(group)
        OPTS.pc = pc;
        OPTS.ps = ps;
        OPTS.B = B;
        OPTS.D = D;
        OPTS.C = C;
        OPTS.sigma = sigma;
        OPTS.lambda = lambda;
        OPTS.weights = weights;
        OPTS.mu = mu;
        OPTS.mueff = mueff;
        OPTS.cc = cc;
       
        OPTS.cs = cs;
        OPTS.c1 = c1;
        OPTS.cmu = cmu;
        OPTS.damps = damps;
        OPTS.chiN = chiN;
        OPTS.countval = countval;
        OPTS.pop = arx';
        group.xmean    = xmean;
        group.bestmem  = arx(:,1)';
        group.bestval  = arfitness(1);
        group.OPTS     = OPTS;
       
        group.delta    = bestval - old_bestval;
        group.cc       = std(arfitness);
        group.iters    = iters;
    end
end