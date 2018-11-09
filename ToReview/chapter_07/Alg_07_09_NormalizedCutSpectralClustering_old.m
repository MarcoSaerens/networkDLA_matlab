function U = Alg_07_09_NormalizedCutSpectralClustering_old(m,A,option)
%   m-way normalized cut spectral clustering algorithm and Ng et al.
%   spectral clustering algorithm
%       Chapter 7 : Clustering nodes
%           Section 8 : Spectral clustering
%               Subsection 4 : Partitioning nodes into three or more
%                              clusters
%
% Arguments:
%
% - m: the number of desired clusters
%
% - A: the n x n adjacency matrix associated to weighted undirected graph G
%
% - option: if normalizedCut, perform a normalized cut spectral clustering
%           if NgNormalisation, perform a normalization as suggested by Ng et al.
%
% Returns:
%
% - U: the n x m binary membership indicator matrix
%
% (c) Simon Gilissen & Edouard Lardinois

% Setting the normalization parameter:

if nargin == 2
    %perform a normalized cut spectral clustering
    normalization = 1;
elseif nargin == 3
    if strcmp(option,'normalizedCut')
        %perform a normalized cut spectral clustering
        normalization = 1;
    elseif strcmp(option,'NgNormalisation')
        %perform a normalization as suggested by Ng et al.
        normalization = 0;
    else
        error('Invalid option name.');
    end
else
    error('Number of arguments is invalid');
end

n = size(A,1);

% Compute the diagonal degree matrix
D = diag(A*ones(n,1));

% Compute unnormalized Laplacian matrix
L = eye(n) - D^(-0.5)*A*D^(-0.5);

% Compute the n x m matrix of eigenvectors of L corresponding to the 
% m smallest non-trivial eigenvalues, sorted by increasing eigenvalues
[V_L , D_L] = eig(L);
[D_L_sorted, idx_D_L_sorted] = sort(D_L, 'ascend');
V_L_sorted = V_L(:,idx_D_L_sorted);
offset = 0;
while D_L_sorted(offset+1) == 0
    offset = offset+1;
end
V = V_L_sorted(:,1+offset:m+offset);
clear offset;

if(normalization==1)
    %loop on nodes
    V = D^(-0.5)*V;
    X = zeros(size(V));
    for i=1:n
        X(i,:) = V(i,:);
    end
    
elseif(normalization==0)
    %loop on nodes
    X = zeros(size(V));
    for i=1:n
        X(i,:) = V(i,:)/ norm(V(i,:))^2 ; %normalizing the vectors
                                          %as suggested by Ng et al.
    end
    
end

% Apply k-means clusters algorithm
y = pdist(X);
Delta = squareform(y);
U = Alg_07_01_StandardkMeansClustering(Delta, m);

end


