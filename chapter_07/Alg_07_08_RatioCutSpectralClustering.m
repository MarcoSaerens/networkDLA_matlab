function U = Alg_07_08_RatioCutSpectralClustering(A, m)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Simon Gilissen & Edouard Lardinois, 
    %          revised by Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description:  Compute the m-way ratio cut spectral clustering 
    %               with optional normalization.
    %
    % INPUT:
    % -------
    % - A : the (n x n) adjacency matrix of a weighted, undirected graph
    % - m : the desired number of clusters
    %
    % OUTPUT:
    % -------
    % - U : the (n x 2) cluster membership matrix with u_ik = 1 if node i 
    %       belongs to cluster k, zero otherwise
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
    
    % The vectors of ones
    e = ones(n, 1);

    % Compute the diagonal degree matrix
    Diag_d = diag(A*e);

    % Compute unnormalized Laplacian matrix
    L = Diag_d - A;

    % Compute the n x m matrix of eigenvectors of L corresponding to the 
    % m smallest non-trivial eigenvalues, sorted by increasing eigenvalues
    [V_L , D_L] = eig(L);
    [~, index_eig] = sort(D_L, 'ascend');
    V_L = V_L(:, index_eig);
    V_L = V_L(:, 1:m);
    
    X = zeros(n, m);
    for i = 1:n
        x_i = V_L(i, :);
        norm_x_i = sqrt(sum(x_i.^2));
        X(i, :) = x_i / norm_x_i;
    end

    % Apply k-means clusters algorithm
    D = zeros(n);
    for i = 1:n
        for j = 1:n
            D(i, j) = (X(i, :) - X(j, :)) * (X(i, :) - X(j, :))';
        end
    end
    U = Alg_07_01_StandardkMeansClustering(D, m);
end
