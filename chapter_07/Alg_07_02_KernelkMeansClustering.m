function U = Alg_07_02_KernelkMeansClustering(K, m)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the kernel k-means clustering
    %
    % INPUT:
    % -------
    % - K : the (n x n) symmetric positive semidefinite kernel matrix
    % - m : a integer containing the number of clusters
    %
    % OUTPUT:
    % -------
    % - U : the (n x m) cluster membership matrix with u_ik = 1 if node i 
    %       belongs to cluster k, zero otherwise
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Checks of arguments
    
    % Check if K is a squared, symmetric matrix
    if ~isequal(K, K')
        error('The kernel matrix is not squared or not symmetric.');
    end
    
    % Check if m is greater than 0
    if m <= 0
        error('The number of groups must be greater than 0');
    end
    
    % Check if m is integer
    if floor(m) ~= m
        error('The number of groups must be an integer');
    end
    
    %% Algorithm 
    
    % Number of nodes, identity matrix, vectors of ones for groups
    [n, ~] = size(K);
    I = eye(n);
    e_m = ones(m, 1);

    % Initialization of a randomly sample m prototypes among the nodes.
    % and creates H matrix containing h_k vectors on columns
    q = randperm(n, m);
    H = zeros(n, m);
    for k = 1:m
        H(q(k), k) = 1;
    end

    % Initialization of the cluster allocation
    U_prev = zeros(n, m);
    U = ones(n, m);
    
    % Main loop
    while ~isequal(U, U_prev)
        
        % Store last cluster allocation
        U_prev = U;
        
        % Compute new cluster allocation
        U = zeros(n, m);
        for i = 1:n
            e_i = I(:, i);
            H_k_e_i = H - e_i * e_m';
            [~, k_star] =  min(diag(H_k_e_i' * K * H_k_e_i));
            U(i, k_star) = 1;
        end

        % Compute prototypes
        n_clust = sum(U, 1);
        H = U ./ n_clust;
        
    end
    
end