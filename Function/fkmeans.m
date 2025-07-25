function [label, centroid, dis] = fkmeans(X, k, options)
% 返回簇标签、簇中心、总平方距离

n = size(X,1);

% option defaults
weight = 0; % uniform unit weighting
careful = 0;% random initialization

if nargin == 3
    if isfield(options, 'weight')
        weight = options.weight;
    end
    if isfield(options,'careful')
        careful = options.careful;
    end
end

% If initial centroids not supplied, choose them
if isscalar(k)
    % centroids not specified
    if careful
        k = spreadseeds(X, k);
    else
        k = X(randsample(size(X,1),k),:);
    end
end

% generate initial labeling of points
[~,label] = max(bsxfun(@minus,k*X',0.5*sum(k.^2,2)));
k = size(k,1); % 实际簇的数量

last = 0;

if ~weight
    % code defactoring for speed
    while any(label ~= last)
        % remove empty clusters
        [~,~,label] = unique(label);
        % transform label into indicator matrix
        ind = sparse(label,1:n,1,k,n,n);
        % compute centroid of each cluster
        centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
        % compute distance of every point to each centroid
        distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
        % assign points to their nearest centroid
        last = label;
        [~,label] = max(distances);
        label = label';
    end
    dis = ind*(sum(X.^2,2) - 2*max(distances)');
else
    while any(label ~= last)
        % remove empty clusters
        [~,~,label] = unique(label);
        % transform label into indicator matrix
        ind = sparse(label,1:n,weight,k,n,n);
        % compute centroid of each cluster
        centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
        % compute distance of every point to each centroid
        distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
        % assign points to their nearest centroid
        last = label;
        [~,label] = max(distances);
        label = label';
    end
    dis = ind*(sum(X.^2,2) - 2*max(distances)');
end
label = label';

function D = sqrdistance(A, B) % 计算两个数据矩阵之间平方欧氏距离
% Square Euclidean distances between all sample pairs
% A:  n1 x d data matrix
% B:  n2 x d data matrix
% WB: n2 x 1 weights for matrix B
% D: n2 x n1 pairwise square distance matrix
%    D(i,j) is the squared distance between A(i,:) and B(j,:)
% Written by Michael Chen (sth4nth@gmail.com). July 2009.
n1 = size(A,1); n2 = size(B,2);
m = (sum(A,1)+sum(B,1))/(n1+n2);
A = bsxfun(@minus,A,m);
B = bsxfun(@minus,B,m);
D = full((-2)*(A*B'));
D = bsxfun(@plus,D,full(sum(B.^2,2))');
D = bsxfun(@plus,D,full(sum(A.^2,2)))';
end

function [S, idx] = spreadseeds(X, k) % 生成初始聚类中心
% X: n x d data matrix
% k: number of seeds
% reference: k-means++: the advantages of careful seeding.
% by David Arthur and Sergei Vassilvitskii
% Adapted from softseeds written by Mo Chen (mochen@ie.cuhk.edu.hk), 
% March 2009.
[n,d] = size(X);
idx = zeros(k,1);
S = zeros(k,d);
D = inf(n,1);
idx(1) = ceil(n.*rand);
S(1,:) = X(idx(1),:);
for i = 2:k
    D = min(D,sqrdistance(S(i-1,:),X));
    idx(i) = find(cumsum(D)/sum(D)>rand,1);
    S(i,:) = X(idx(i),:);
end
end

end
