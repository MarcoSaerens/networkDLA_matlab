function U = Alg_07_11_LatentClassClustering(A, m)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description:  Compute the latent class clustering of a graph
    %
    % INPUT:
    % -------
    % - A : the (n x n) adjacency matrix of a weighted, undirected graph
    % - m : the desired number of clusters
    %
    % OUTPUT:
    % -------
    % - U : the (n x m)  membership matrix
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Checks of arguments
    
    [n, ~] = size(A);
    
    % Check if A is a squared, symmetric matrix
    if ~isequal(A, A')
        error('The adjacency matrix is not squared or not symmetric.');
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
    
    % The vector of ones
    e = ones(n, 1);
    
    % Prob array 
    Prob_ij_k = (ones(n, n, m) - 1e-2 * rand(n, n, m) ) / m ;
    for i = 1:n
        for j = 1:n
            Prob_ij_k(i, j, :) =  Prob_ij_k(i, j, :) ./ sum( Prob_ij_k(i, j, :) );
        end
    end
    
    Prob_ij_k_old = Prob_ij_k + 100;
    
    p_k = zeros(n);
    p_i_k = zeros(n, m);
    while ~isequal(round(Prob_ij_k_old, 4), round(Prob_ij_k, 4))
        
        % Store last value
        Prob_ij_k_old = Prob_ij_k;
        
        for k = 1:m 
            
            %%%% M STEP %%%% 
            p_k(k) = sum(sum(Prob_ij_k(:, :, k) .* A)) / sum(sum(A));
            p_i_k(:, k) = sum(Prob_ij_k(:, :, k) .* A, 2);
            p_i_k(:, k) = p_i_k(:, k) / sum(p_i_k(:, k));
            p_j_k = sum(Prob_ij_k(:, :, k) .* A, 1);
            p_j_k = p_j_k / sum(p_j_k);
            
            %%%% E STEP %%%% 
            Prob_ij_k(:, :, k) = p_k(k) .* (p_i_k(:, k) * e') .* (e * p_j_k); 
        end
        
        % Normalization 
        for i = 1:n
            for j = 1:n
                Prob_ij_k(i, j, :) =  Prob_ij_k(i, j, :) ./ sum( Prob_ij_k(i, j, :) );
            end
        end
    end
    
    U = zeros(n, m);
    for i = 1:m
        for k = 1:m
            U(i, k) = p_i_k(i, k) * p_k(k);
        end
        U(i, :) = U(i, :) / sum(U(i, :));
    end
    
end
