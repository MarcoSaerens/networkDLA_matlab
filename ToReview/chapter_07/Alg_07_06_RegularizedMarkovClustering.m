function U = Alg_07_06_RegularizedMarkovClustering(A, r)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the regularized Markov clustering of nodes
    %
    % INPUT:
    % -------
    % - A : the (n x n) adjacency matrix of a weighted, undirected graph
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
    % Initial matrix 
    Q = P;
    
    % Main loop
    Q_prev = Q + 1000;
    while ~isequal(Q, Q_prev)
        % Store last value
        Q_prev = Q;
        % Expansion
        Q = P*Q;
        % Inflation
        Q = Q.^r;
        % Normalization
        Q = (diag(Q*e))^(-1) * Q;
    end
    
    % Pruning
    U = Q(:, sum(Q, 1) > 0);
end

