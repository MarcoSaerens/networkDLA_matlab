function X = Alg_10_05_LatentSocialMap(A, p, X_0)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Authors: Gilissen & Edouard Lardinois, revised by Guillaume Guex (2017).
    %
    % Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
    %         "Algorithms and models for network data and link analysis". 
    %         Cambridge University Press.
    %
    % Description: Computes the latent social embedding of a graph.
    %
    % INPUT:
    % -------
    % - A : a (n x n) unweighted adjacency matrix associated with an undirected
    %       graph G.
    % - p : a integer containing the number of dimensions kept for 
    %       the embedding.
    % - X_0 : a (n x p_0) initial position for the n nodes.
    %
    % OUTPUT:
    % -------
    % - X : the (n x p) data matrix containing the coordinates of the nodes 
    %       for the embedding.    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Checks of arguments

    % Check if squared matrix
    [n, m] = size(A);
    if n ~= m
        error('The adjacency matrix must be a square matrix')
    end

    % Check if symmetric matrix / graph is undirected
    if ~isequal(A, A')
        error('The adjacency matrix is not symmetric.')
    end

    % Check if binary matrix
    if ~isequal(A, (A > 0))
        error('The adjacency matrix is not binary.')
    end

    % Check the dimension of Initial matrix
    [n_0, p_0] = size(X_0);
    if n_0 ~= n
        error('The initial position matrix must have a number of rows equal to the size of the adjacency matrix')
    end

    % Check if number of dimensions is pertinent
    if p > p_0
        error('The number of kept dimensions must be less or equal than to the number of dimensions');
    end

    %% Algorithm 

    % Setting initial parameters
    alpha = 1;
    X = X_0(:,1:p);
    D = zeros(n);
    Y_hat = zeros(n);
    B = eye(n*p + 1);
    gradient_unfold_prev = false;
    minus_l_prev = realmax;
    minus_l = 1e15;
    
    % Main loop
    while abs(minus_l_prev - minus_l) > 1e-4
        
        % Store previous minus likelihood
        minus_l_prev = minus_l;
        
        for i = 1:n
            for j = 1:n
                % Compute the distances in the latent space
                D(i, j) = sqrt( (X(i, :) - X(j, :)) * (X(i, :) - X(j, :))' ); 
                % Compute the predicted pobability of each link
                Y_hat(i, j) = 1 / ( 1 + exp( -alpha * (1 - D(i, j)) ) ); 
            end
        end

        % Compute the gradient with respect to alpha    
        G_alpha_mat = (D - 1) .* (Y_hat - A);
        G_alpha_mat(1:n+1:end) = 0;
        gradient_alpha = 1/2 * sum(sum(G_alpha_mat));

        % Compute the gradient with respect to position vector x_k
        Gradient_X = zeros(n, p);
        for k = 1:n
            for j = 1:n
                if j ~= k
                    xk_xj = (X(k, :) - X(j, :));
                    Gradient_X(k, :) = Gradient_X(k, :) + alpha * (Y_hat(k, j) - A(k, j)) * xk_xj / (norm(xk_xj) + 1e-20) ;
                end
            end
        end
        
        %%% Perform a min search of a BFGS quasi-Newton step on minus likelihood
               
        % Unfold gradient
        gradient_unfold = [gradient_alpha; Gradient_X(:)];
        
        % Compute B (if this is not the first step)
        if gradient_unfold_prev ~= false
            y = gradient_unfold - gradient_unfold_prev;
            B = B + y*y' / (y'*dir_unfold) - (B*dir_unfold)*(dir_unfold'*B) / (dir_unfold'*B*dir_unfold);
        end
        
        % Obtain direction
        dir_unfold = linsolve(B, -gradient_unfold);
        
        % Fold
        dir_alpha = dir_unfold(1);
        Dir_X = reshape(dir_unfold(2:end), [n, p]);
        
        % Find stepsize
        [beta_opt, ~] = fminsearch(@(beta) minus_likelihood(A, alpha, X, dir_alpha, Dir_X, beta), 0.5);
        
        % Make step
        X = X + beta_opt*Dir_X;
        alpha = alpha + beta_opt*dir_alpha;
        
        % Store new minus likelihood
        minus_l = minus_likelihood(A, alpha, X, dir_alpha, Dir_X, 0);
    end
end

%% Minus likelihood function

function minus_l = minus_likelihood(A, alpha, X, dir_alpha, Dir_X, beta)
    n = size(X, 1);
    
    X_beta = X + beta*Dir_X;
    alpha_beta = alpha + beta*dir_alpha;

    D = zeros(n);
    for i=1:n
        for j=1:n
            D(i,j) = sqrt( (X_beta(i,:)-X_beta(j,:))*(X_beta(i,:)-X_beta(j,:))' );
        end
    end
    
    components_of_likelihood = (alpha_beta * A .* (1 - D) - log(1 + exp(alpha_beta * (1 - D))) );
    components_of_likelihood(1:n+1:end) = 0;
    minus_l = - 1/2 * sum(sum(components_of_likelihood));
end
