function K_BoP = Alg_03_08_BagOfPathsCovariance(C, P_ref, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the bag-of-paths covariance matrix between nodes 
%
% INPUT:
% -------
% - C : the (n x n) adjacency matrix associated to a weighted, strongly
%       connected graph,containing affinities.
% - P_ref : the (n x n) reference transition probabilities matrix.
% - theta : a strictly positive parameter controlling the degree of 
%           randomness.
%
% OUTPUT:
% -------
% - K_BoP : the bag of paths covariance matrix between pairs of nodes
%       containing the elements K_BoP(k, l) = cov(k, l)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments 

% Check if C is squared  
[n, m] = size(C);
if n ~= m 
    error('The cost matrix is not squared')
end

% Check if P_ref has the same dimension as C 
[o, p] = size(P_ref);
if (o ~= n) || (p ~= m)
    error('The reference transition probabilities matrix does not correspond to the cost matrix')
end

% Check if theta is positive
if theta <= 0
    error('The randomness parameter is not stricly positive')
end

%% Algorithm 

% Elemental exponential and multiplication 
W = P_ref .* exp(-theta*C);

% The fundamental matrix 
Z = (eye(n) - W)^(-1);

% Sums of cols, rows, and all elements  
z_dotk = sum(Z, 1);
z_kdot = sum(Z, 2);
z_dotdot = sum(z_dotk);

% The partition function 
z_bar = z_dotdot - n;

K_BoP = zeros(n, n);
for k = 1:n
    for l = k:n
        if l == k
            delta = 1;
        else
            delta = 0;
        end
        p_1 = (z_dotk(k) - 1) * z_kdot(k) * delta;
        p_2 = z_dotk(k) * (z_kdot(l) - 1) * (Z(l, k) - delta);
        p_3 = z_kdot(l) * (z_dotk(k) - 1) * (Z(k, l) - delta);
        p_4 = (1/z_bar) * ( z_kdot(k) * z_kdot(l) * (z_dotk(k)-1) * (z_dotk(l)-1) );
        
        K_BoP(k, l) = (1/z_bar) * (p_1 + p_2 + p_3 - p_4);

        
        K_BoP(l, k) = K_BoP(k, l);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%