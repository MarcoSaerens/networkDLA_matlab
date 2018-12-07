function U = Alg_08_07_LouvainMethod(A,mix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Ilkka Kivimaki (2015),revised by Sylvain Courtain (2017).
% Direct source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
% "Algorithms and models for network data and link analysis".
% Cambridge University Press.
%
% Performs the Louvain Method
%
% INPUT:
% ------
% - A weighted undirected graph containing n nodes
% - A, the n x n adjancencty matrix of the graph
%
% OUTPUT:
% -------
% - U, the n x m membership matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECK ERRORS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    error('The algorithm needs an adjacency matrix as input!');
end

[N,n] = size(A);
if (N ~= n)
    error('Non-square matrix A!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1 : Local optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy the matrix A before mix
Adj = A;
% Mix the nodes in random order:
if mix
    N = size(A, 1);
    randidx = randperm(N);
    A = A(randidx, :);
    A = A(:, randidx);
else
    randidx = 1:N;
end
level=1;
% Initialize the binary membership matrix
N = size(A, 1);  % N = # of nodes
assigns = 1:N;    % Assign each node initially to its own cluster.
U = eye(N);   % Cluster membership matrix


% Compute the modularity matrix
Q = ComputeModularityMatrix(A);

i=0;      % Initialize node index
conv=1;   % Initialize convergence counter

while conv<N   % Repeat iterations until no move of a node occurs
    
    i = i+1;    % Choose node i
    clust_i = assigns(i);   % Which cluster does i belong to?
    U_clust_i = U(:, clust_i);   % Membership vector of i's cluster
    
    % Contribution to modularity of having i in clust_i:
    Qloss = U_clust_i'*Q(:, i) + Q(i, :)*U_clust_i - 2*Q(i, i);
    
    % Compute the set of all ajacent clusters of node i
    neigh_i = (Adj(i, :) > 0);   % Neighbors of i (according to A)
    clust_neigh_i = unique(assigns(neigh_i));    % Cluster nhood of i
    clust_neigh_i = setdiff(clust_neigh_i, clust_i); % Don't consider the cluster of i
    
    % Combine i with all neighboring clusters and find the best:
    best_Qgain = -Inf;   % Initialize best modularity change
    
    for J=clust_neigh_i
        
        % Contribution to modularity of having i in cluster J:
        Qgain = U(:,J)'*Q(:, i) + Q(i, :)*U(:, J);
        
        % Find the best move for node i
        if Qgain > best_Qgain    % Check if better than earlier
            best_Qgain = Qgain;
            best_J = J;   % Save index
        end
    end
    
    
    % Move i to cluster best_J, if it increases modularity:
    if best_Qgain > max(0, Qloss)
        assigns(i) = best_J;
        U(i,clust_i) = 0;
        U(i,best_J) = 1;
        
        conv = 1;   % Start the countdown from the beginning
        
        % Otherwise increase the counter:
    else
        conv = conv+1;
    end
    
    % After node N, go back to node 1:
    if i==N
        i=0;
    end
end
% Reorganize clusters, rebuild U and calculate cluster sizes:
[U, assigns] = LMBOP_reorganize(assigns);
% If the indexes were mixed, then reorganize them:
if mix
    origassigns = 0*assigns;  % allocation
    for i = 1:N
        origassigns(i) = assigns(randidx==i);
    end
end
fullassigns = assigns';



while true
    %% Step 2 and Convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Construct the new graph by merging clusters:
    A_merge = U'*A*U;
    %% Step 1 : Local optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialize the binary membership matrix
    M = size(A_merge, 1);  % N = # of nodes
    newassigns = 1:M;    % Assign each node initially to its own cluster.
    U = eye(M);   % Cluster membership matrix
    level = level + 1;
    
    % Compute the modularity matrix
    Q = ComputeModularityMatrix(A_merge);
    
    
    
    i=0;      % Initialize node index
    conv=1;   % Initialize convergence counter
    
    while conv<M   % Repeat iterations until no move of a node occurs
        
        i = i+1;    % Choose node i
        clust_i = newassigns(i);   % Which cluster does i belong to?
        U_clust_i = U(:, clust_i);   % Membership vector of i's cluster
        
        % Contribution to modularity of having i in clust_i:
        Qloss = U_clust_i'*Q(:, i) + Q(i, :)*U_clust_i - 2*Q(i, i);
        
        % Compute the set of all ajacent clusters of node i
        neigh_i = (A_merge(i, :) > 0);   % Neighbors of i (according to A)
        clust_neigh_i = unique(newassigns(neigh_i));    % Cluster nhood of i
        clust_neigh_i = setdiff(clust_neigh_i, clust_i); % Don't consider the cluster of i
        
        % Combine i with all neighboring clusters and find the best:
        best_Qgain = -Inf;   % Initialize best modularity change
        
        for J=clust_neigh_i
            
            % Contribution to modularity of having i in cluster J:
            Qgain = U(:, J)'*Q(:, i) + Q(i, :)*U(:, J);
            
            % Find the best move for node i
            if Qgain > best_Qgain    % Check if better than earlier
                best_Qgain = Qgain;
                best_J = J;   % Save index
            end
        end
        
        
        % Move i to cluster best_J, if it increases modularity:
        if best_Qgain > max(0, Qloss)
            newassigns(i) = best_J;
            U(i, clust_i) = 0;
            U(i, best_J) = 1;
            
            conv = 1;   % Start the countdown from the beginning
            
            % Otherwise increase the counter:
        else
            conv = conv+1;
        end
        
        % After node N, go back to node 1:
        if i==M
            i=0;
        end
    end
    [U, newassigns] = LMBOP_reorganize(newassigns);
    %% Convergence Step
    % Announce the clustering in terms of the original graph:
    oldassigns = fullassigns;
    for i = 1:max(oldassigns)
        J = (oldassigns == i);
        fullassigns(J) = newassigns(i);
    end
    % Check convergence:
    if fullassigns == oldassigns
        U = LMBOP_reorganize(origassigns);
        break;
    end
    
    % If not converged, continue:
    
    % Reorganize clusters, rebuild U and compute cluster sizes:
    [U, fullassigns] = LMBOP_reorganize(fullassigns);
    % Reorganize according to original indices:
    origassigns = fullassigns;
    for i = 1:N
        origassigns(i) = fullassigns(randidx==i);
    end
end


function Q = ComputeModularityMatrix(A)
% Compute the modularity matrix
vol = sum(sum(A));
din  = sum(A,1); % the indegree vector
dout = sum(A,2); % the outdegree vector
D  = (dout * din)/vol;
Q  = (A - D); % Compute the modularity matrix

function [U, newassigns] = LMBOP_reorganize(oldassigns)

N = numel(oldassigns);

U = zeros(N);
for i = 1:N,
    U(:,i) = (oldassigns == i);
end;

n_of_clusters = numel(unique(oldassigns));
[clsizes, idx] = sort(sum(U), 'descend');
idx = idx(1:n_of_clusters);
U = U(:,idx);
newassigns = 0*oldassigns;
for i = 1:n_of_clusters
    old_id = idx(i);
    newassigns(oldassigns == old_id) = i;
end;
