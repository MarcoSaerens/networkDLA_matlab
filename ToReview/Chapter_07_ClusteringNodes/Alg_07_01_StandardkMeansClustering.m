function U = Alg_07_01_StandardkMeansClustering(D, m)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the standard distance-based k-means clustering
    %
    % INPUT:
    % -------
    % - D : the (n x n) symmetric distance matrix containing distances
    %       between nodes
    % - m : a integer containing the number of clusters
    %
    % OUTPUT:
    % -------
    % - U : the (n x m) cluster membership matrix with u_ik = 1 if node i 
    %       belongs to cluster k, zero otherwise
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Checks of arguments
    
    % Check if D is a squared, symmetric matrix
    if ~isequal(D, D')
        error('The distance matrix is not squared or not symmetric.');
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
    
    % Number of nodes
    [n, ~] = size(D);

    % Initialization of a randomly sample m prototypes among the nodes.
    q = randperm(n, m);

    % Initialization of values of the objective function
    j_prev = realmax;
    j = 1e15;
    
    % Main loop
    while abs(j - j_prev) > 1e-4
        
        % Store last value of objective function
        j_prev = j;
        
        % Reallocation of nodes
        [~, l] = min(D(q, :));

        % Recomputation of the prototype of each cluster and
        % the objective function
        j = 0;
        for k = 1:m           
            cluster_index = (l==k);
            [min_val, q(k)] = min(sum(D(cluster_index, :), 1));
            j = j + min_val;
        end
    end

    % Fill in cluster membership matrix
    U = zeros(n, m);
    for i = 1:n
        U(i,l(i)) = 1;
    end
    
end