% 基于密度排序
function rho_mean = kernel (temp_x,temp_pop,dim,imp_total)
%     R = 0.01;
   
    h = 1.5;
%     if imp_total ==0
%         h = 2;
%     else
%         h = min(2,1/(1+exp(((imp_total-4))))/0.05);
%     end
%     h = 0.0000001;
%     h = 0.000000000000000000000000000000000001;
    % 阈值
    rho_min = exp(-h^2 / (2*h^2));
 
    if isempty(temp_pop)
        rho_mean = zeros(size(temp_x,1),1);
    else
        dist = pdist2(temp_x,temp_pop);

        rho = exp(-dist.^2 ./ (2 * h^2));
        rho(rho < rho_min) = 0;
        rho_mean = mean(rho, 2);
    end
end