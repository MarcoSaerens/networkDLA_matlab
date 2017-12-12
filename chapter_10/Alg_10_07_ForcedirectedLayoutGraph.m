function X = Alg_10_07_ForcedirectedLayoutGraph(W, w, X_0, a, r)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes (a, r) force-directed layout for a graph
    %
    % INPUT:
    % -------
    % - W:  a n x n matrix containing edges weights a weighted, 
    %       undirected graph G.
    % - w:  a n x 1 vector containing non-negative node weights
    % - X_0: n x p initial position for the n nodes
    % - a:  the attraction power parameter
    % - r:  the repulsion power parameter
    %
    % OUTPUT:
    % -------
    % - X : the (n x p) data matrix containing the coordinates of the nodes 
    %       on its rows.    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Checks of arguments
    
    [n, ~] = size(W);
    
    % Check if W is a symmetric matrix
    if ~isequal(W, W')
        error('The edge weights matrix is not symmetric.')
    end
    
    % Check if W are positive
    if min(min(W)) < 0
        error('The edge weights matrix must be a non-negative.')
    end
    
    % Check if w are positive
    if min(min(w)) < 0
        error('The node weights vector must be a non-negative.')
    end
    
    % Check if w correspond to W
    n_w = length(w);
    if n ~= n_w
        error('The node weight vector length must correspond to the edge weights matrix')
    end

    % Check the dimension of Initial matrix
    [n_0, p] = size(X_0);
    if n ~= n_0
        error('Input X_0 must have a number of rows equal to the size of W')
    end
    
    %% Algorithm 
    
    % Initialize the embedding
    X = X_0;
    energy_prev = realmax;
    energy = 1e15;
    B = eye(n*p);
    gradient_unfold_prev = false;
    
    while abs(energy_prev - energy) > 1e-6
        
        % Store previous energy value
        energy_prev = energy;
        
        % Compute the gradient
        Gradient_X = zeros(n, p); 
        random_index = randperm(n);
        for k = random_index
            for j = 1:n
                if j ~= k
                    xk_xj = X(k, :) - X(j, :);
                    norm_xk_xj = (norm(xk_xj) + 1e-20);
                    Gradient_X(k, :) = Gradient_X(k, :) + W(j, k)*(norm_xk_xj^(a-1) - w(j)*w(k)*norm_xk_xj^(r-1))*xk_xj;
                end
            end
        end
        
        %%% Perform a min search of a BFGS quasi-Newton step on minus likelihood
               
        % Unfold gradient
        gradient_unfold = Gradient_X(:);
        
        % Compute B (if this is not the first step)
        if gradient_unfold_prev ~= false
            y = gradient_unfold - gradient_unfold_prev;
            B = B + y*y' / (y'*dir_unfold) - (B*dir_unfold)*(dir_unfold'*B) / (dir_unfold'*B*dir_unfold);
        end
        
        % Obtain direction
        dir_unfold = linsolve(B, -gradient_unfold);
        
        % Fold
        Dir_X = reshape(dir_unfold, [n, p]);
        
        % Find stepsize
        [beta_opt, ~] = fminsearch(@(beta) potential_energy(W, w, X, Dir_X, a, r, beta), 1e-2);
        
        % Make step
        X = X + beta_opt*Dir_X;
        
        % Store new energy
        energy = potential_energy(W, w, X, Dir_X, a, r, 0);
        
    end
end

function e = potential_energy(W, w, X, Dir_X, a, r, beta)
    [n, ~] = size(X);
    X_beta = X + beta*Dir_X;
    e = 0;
    for i = 1:(n-1)
        for j = (i+1):n
            norm_xi_xj = norm(X_beta(i, :) - X_beta(j, :)) + 1e-20;
            if a == -1
                e = e + W(i, j) * log(norm_xi_xj);
            else
                e = e + W(i, j) * norm_xi_xj^(a + 1) / (a + 1);
            end
            if r == -1
                e = e - w(i)*w(j) * log(norm_xi_xj);
            else
                e = e - w(i)*w(j) * norm_xi_xj^(r + 1) / (r + 1);
            end
        end
    end
end





