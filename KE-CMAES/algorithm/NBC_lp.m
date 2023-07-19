function [species, meandis]= NBC_lp(pop)
% function [nbc]= NBC_large_population(pop)
%NBC 此处显示有关此函数的摘要
%   此处显示详细说明
[n, dim] = size(pop);       % 边
min_num_edge = 10;

factor = 4 - log(dim);

nbc=zeros(n,3);             % 信息矩阵，起始点，终止点，距离
nbc(1:n, 1) = 1:n;          % 起始点
nbc(1, 2) = -1;             % 终止点
nbc(1, 3) = 0;              % 距离

for i = 2 : n
    arrdis = pdist2(pop(i, :), pop(1:i-1, :));
    [u, v] = min(arrdis);
    nbc(i,2) = v;
    nbc(i,3) = u;
end

% filename = sprintf('./nbc/nbc_F%d_runs%d.mat', func, runs);
% load(filename);

% nbc = nbc(1:n, :);

% meandis=factor*mean(nbc(2:n,3));
% nbc(nbc(:,3)>meandis,2)=-1;
% nbc(nbc(:,3)>meandis,3)=0;
meandis=factor*mean(nbc(2:n,3));
if  length(nbc(nbc(:,3)>meandis)) >= min_num_edge
    nbc(nbc(:,3)>meandis,2)=-1;
    nbc(nbc(:,3)>meandis,3)=0;
else
    pred = sortrows(nbc,3,'descend');
    meandis = pred(min_num_edge, 3);
    nbc(nbc(:,3)>=meandis,2)=-1;
    nbc(nbc(:,3)>=meandis,3)=0;
end


seeds=nbc(nbc(:,2)==-1,1);
m=zeros(n,2); % 保存各个粒子索引以及对应的簇号
m(1:n,1)=1:n;
for i=1:n
    j=nbc(i,2);
    k=j;
    while j~=-1
        k=j;
        j=nbc(j,2);
    end
    if k==-1
        m(i,2)=i;
    else
        m(i,2)=k;
    end
end % 此时每个粒子对应的是该粒子所在簇的seed编号

% construct the result
    species = struct();
    for i=1:length(seeds)
       species(i).seed = seeds(i);
       species(i).idx = m(m(:, 2) == seeds(i), 1);
       species(i).len = length(species(i).idx);
    end

end

