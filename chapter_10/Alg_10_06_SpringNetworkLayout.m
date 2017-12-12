function X = Alg_10_06_SpringNetworkLayout(D, X_0, l_0, k)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the spring network layout of a graph.
    %
    % INPUT:
    % -------
    % - D:  n x n symmetric all-pairs shortest-path distances matrix 
    %       associated with graph G, providing distances Delta_ij 
    %       for each pair of nodes i, j
    % - X_0:    n x p initial position for the n nodes
    % - l_0: drawing length constant
    % - k:  global strength constant
    %
    % OUTPUT:
    % -------
    % - X : the (n x p) data matrix containing the coordinates of the nodes 
    %       on its rows.    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Checks of arguments
    
    [n, ~] = size(D);

    % Check if D is a symmetric matrix
    if ~isequal(D, D')
        error('The adjacency matrix is not symmetric.')
    end

    % Check the dimension of Initial matrix
    [n_0, p] = size(X_0);
    if n ~= n_0
        error('Input X_0 must have a number of rows equal to the size of D')
    end

    % Check if display length is positive
    if l_0 <= 0
        error('Display length l_0 must be more than zero')
    end

    % Check if global strength constant is positive
    if k <= 0
        error('Global strength constant k must be more than zero')
    end
    
    %% Algorithm 

    % Set equilibrium length of each spring
    L_0 = l_0 * D / max(max(D));

    % Set stiffness of each spring
    K_stiff = k ./ D.^2 ;

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
        for k = 1:n
            for j = 1:n
                if j ~= k
                    xk_xj = X(k, :) - X(j, :);
                    Gradient_X(k, :) = Gradient_X(k, :) + K_stiff(k,j) * (1 - L_0(k,j)/(norm(xk_xj) + 1e-20)) * xk_xj;
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
        [beta_opt, ~] = fminsearch(@(beta) energy_spring(K_stiff, L_0, X, Dir_X, beta), 1e-3);
        
        % Make step
        X = X + beta_opt*Dir_X;
        
        % Store new energy
        energy = energy_spring(K_stiff, L_0, X, Gradient_X, 0);
    end
end

%% Energy spring function

function e = energy_spring(K_stiff, L_0, X, Dir_X, beta)
    [n, ~] = size(X);
    X_beta = X + beta*Dir_X;
    e = 0;
    for i = 1:(n-1)
        for j = (i+1):n
            e = e + K_stiff(i, j)/2 * ( norm(X_beta(i, :)-X_beta(j, :)) - L_0(i, j) )^2;
        end
    end
end