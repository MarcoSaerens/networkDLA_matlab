function U = Alg_07_04_LabelPropagationClustering(A, lambda)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the simple label propagation clustering
    %
    % INPUT:
    % -------
    % - A : the (n x n) adjacency matrix of a weighted, undirected graph
    % - lambda: the parameter (>= 0) balancing the size of the clusters 
    %
    % OUTPUT:
    % -------
    % - U : the (n x m) cluster membership matrix with u_ik = 1 if node i 
    %       belongs to cluster k, zero otherwise
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Checks of arguments
    
    % Check if A is a squared, symmetric matrix
    if ~isequal(A, A')
        error('The adjacency matrix is not squared or not symmetric.');
    end
    
    % Check if lambda is greater or equal to 0
    if lambda < 0
        error('The balancing parameter is not non-negative');
    end

    %% Algorithm 
    
    % Number of nodes
    [n, ~] = size(A);
    
    % Initialization of the cluster allocation
    U = eye(n);
    l = 1:n;
    n_clust = ones(n, 1);
    m = n;
    
    % Main loop
    U_prev = zeros(n, m);
    while ~isequal(U, U_prev)
        
        % Store last cluster allocation
        U_prev = U;
        
        for k = 1:n
            
            %%% Try to move node k from current cluster l(k) to a new cluster
            
            % Current value of the term involving node k in H_prime
            h_k_prime = sum( A(k, :) .* (l == l(k)) ) - lambda*( n_clust(l(k)) - 1 );
            
            % Maximum value achieved when changing the cluster label 
            h_k_prime_new = zeros(m, 1);
            for clust = 1:m
                h_k_prime_new(clust) = sum( A(k, :) .* (l == clust) ) - lambda*( n_clust(clust) - (clust == l(k)) );
            end
            [h_k_prime_star, l_star] = max(h_k_prime_new);
            
            % Test if the relabeling increases h
            if (h_k_prime_star > h_k_prime) && (sum(l == l_star) > 0 ) 
                n_clust(l_star) = n_clust(l_star) + 1;
                n_clust(k) = n_clust(k) - 1;
                U(k, l_star) = 1;
                U(k, l(k)) = 0;
                l(k) = l_star;
            end
            
        end 
        
    end
    U = U(:, sum(U, 1) > 0);
    
end