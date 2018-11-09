function U = Alg_07_03_IterativeKernelkMeansClustering(K, m)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the iterative kernel k-means clustering
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

    % Compute cluster allocation
    U = zeros(n, m);
    l = zeros(n, 1);
    for i = 1:n
        e_i = I(:, i);
        H_k_e_i = H - e_i * e_m';
        [~, k_star] =  min(diag(H_k_e_i' * K * H_k_e_i));
        U(i, k_star) = 1;
        l(i) = k_star;
    end

    % Compute prototypes
    n_clust = sum(U, 1);
    H = U ./ n_clust;
    
    % Initialization of the cluster allocation
    U_prev = zeros(n, m);
    
    % Main loop
    while ~isequal(U, U_prev)
        
        % Store last cluster allocation
        U_prev = U;
        
        for i = 1:n
            e_i = I(:, i);
            delta_j_i = zeros(m, 1);
            for k = 1:m
                % Try to move node i from current cluster to cluster k
                h_k_e_i = H(:, k) - e_i;
                h_l_e_i = H(:, l(i)) - e_i;
                delta_j_i(k) = n_clust(k)*(h_k_e_i'*K*h_k_e_i)/(n_clust(k) + 1) - n_clust(l(i))*(h_l_e_i'*K*h_l_e_i)/(n_clust(l(i)) - 1);  
            end
            
            % Lowest delta_j index
            [~, k_star] = min(delta_j_i);
            
            % Test if the move decrease the within cluster inertia
            if delta_j_i(k_star) < 0
                % Remove from l(i)
                H(:, l(i)) = ( n_clust(l(i)) * H(:, l(i)) - e_i ) / (n_clust(l(i)) - 1);
                n_clust(l(i)) = n_clust(l(i)) + 1;
                U(i, l(i)) = 0;
                % Add to k_star
                H(:, k_star) = ( n_clust(k_star) * H(:, k_star) + e_i ) / (n_clust(l(i)) + 1);
                n_clust(k_star) = n_clust(k_star) + 1;
                U(i, k_star) = 1;
                % Update l
                l(i) = k_star;
            end       
        end
        
    end
    
end