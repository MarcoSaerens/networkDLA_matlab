function U = Alg_07_07_KernighanLinClustering(A, U_init)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description:  A heuristic procedure greedily improving graph cut:
    %               Kernighan-Lin algorithm.
    %
    % INPUT:
    % -------
    % - A : the (n x n) adjacency matrix of a weighted, undirected graph
    % - U_init : the initial (n x 2) membership matrix for the 2 groups
    %
    % OUTPUT:
    % -------
    % - U : the (n x 2) cluster membership matrix with u_ik = 1 if node i 
    %       belongs to cluster k, zero otherwise
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Checks of arguments
    
    [n, ~] = size(A);
    
    % Check if A is a squared, symmetric matrix
    if ~isequal(A, A')
        error('The adjacency matrix is not squared or not symmetric.');
    end
    
    % Check if U_init has the same number of rows as A
    [n_u, m_u] = size(U_init);
    if n_u ~= n
        error('The U_init matrix has not the same number of rows as A.');
    end
    
    % Check if U_init has only 2 groups
    if m_u ~= 2 
        error('The U_init matrix must contain 2 column.');
    end

    %% Algorithm 
    
    % The identity matrix, vectors of ones
    I = eye(n);
    e = ones(n, 1);
    
    % The two cluster vectors
    u_1 = U_init(:, 1);
    u_2 = U_init(:, 2);
    
    % Main loop
    convergence = false;
    while ~convergence
        
        % Compute the initial graph cut 
        J_0 = u_1' * A * u_2;
        % Compute the initial dif
        dif = A * (u_1 - u_2);
        
        % Initialize the set of candidate nodes and cluster 
        ind_c1 = u_1;
        ind_c2 = u_2;
        
        % Determine min cluster
        m = min( sum(u_1) , sum(u_2));
          
        % Init vectors
        i_star_vec = zeros(m, 1);
        j_star_vec = zeros(m, 1);
        J_vec = zeros(m + 1, 1);
        J_vec(1) = J_0;
        
        % Try sequentially to swap nodes
        for tau = 1:m
            % Find best candidate swap
            Difi_difj_Aij = (dif .* ind_c1) * e' - e * (dif .* ind_c2)' + A;
            [max_col, i_vec] = max(Difi_difj_Aij);
            [delta_J, j_star_vec(tau)] = max(max_col);
            i_star_vec(tau) = i_vec(j_star_vec(tau));
            
            % Update membership vectors 
            ind_c1(i_star_vec(tau)) = ind_c1(i_star_vec(tau)) - 1; 
            ind_c2(j_star_vec(tau)) = ind_c2(j_star_vec(tau)) - 1; 
            
            % Update graph cut 
            J_vec(tau + 1) = J_vec(tau) + delta_J;
            
            % Update dif
            dif = dif + 2*A*( I(:, j_star_vec(tau)) - I(:, i_star_vec(tau)) );         
        end
        
        % Find the step t corresponding to the best graph cut 
        [J_t, t] = min(J_vec);
        t = t - 1;
        
        % Test if there is an improvement
        if t > 0
            % Update u_1 and u_2
            u_1(i_star_vec(1:t)) = 0;
            u_1(j_star_vec(1:t)) = 1;
            u_2(j_star_vec(1:t)) = 0;
            u_2(i_star_vec(1:t)) = 1;
        end
        
        % Check convergence 
        if J_0 == J_t
            convergence = true;
        end
        
    end
    
    % Concanate
    U = [u_1, u_2];
end

