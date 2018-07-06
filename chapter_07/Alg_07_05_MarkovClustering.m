function U = Alg_07_05_MarkovClustering(A, q, r)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the Markov clustering of nodes
    %
    % INPUT:
    % -------
    % - A : the (n x n) adjacency matrix of a weighted, undirected graph
    % - q : the expension parameter (q > 1)
    % - r : the inflation parameter (r > 1)
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
    
    % Check if q is greater than 1
    if q <= 1 
        error('The expansion parameter q is not greater than 1');
    end
    
    % Make sure q in integer
    q = floor(q);
    
    % Check if r is greater than 1
    if r <= 1 
        error('The inflation parameter r is not greater than 1');
    end

    %% Algorithm 
    
    % Number of nodes, identity matrix, vectors of ones
    [n, ~] = size(A);
    I = eye(n);
    e = ones(n, 1);
    
    % Add small self loop
    A = A + 1e-4 * I;
    % Degree matrix
    Diag_d = diag(A*e);
    % Transition matrix
    P = (Diag_d)^(-1) * A;
    
    % Main loop
    P_prev = P + 1000;
    while ~isequal(P, P_prev)
        % Store last value
        P_prev = P;
        % Expansion
        P = P^q;
        % Inflation
        P = P.^r;
        % Normalization
        P = (diag(P*e))^(-1) * P;
    end
    
    % Pruning
    U = P(:, sum(P, 1) > 0);
end

