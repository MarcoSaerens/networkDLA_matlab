function cr = Alg_04_16_BagOfPathsCriticality(A, C, theta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Bertrand Lebichot, revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016),
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%         See also: Lebichot B. and Saerens M. (2018), "A bag-of-paths
%         network criticality measure". Neurocomputing, 275, pp. 224-236.
%
% Description: The bag of paths approach for computing an approximate
% criticality measure on the nodes of a weighted directed or undirected,
% strongly connected, graph G without self-loops. See the Neurocomputing
% paper for more information.
%
% INPUT:
% ------
% - A: the n x n adjacency matrix of a weighted undirected or directed,
%   strongly connected, graph G containing n nodes.
% - C: the n x n cost matrix C associated to G (if not specified,
%   the costs are the inverse of the affinities, but other 
%   choices are possible).
% - theta: the (scalar) inverse temperature parameter.
%
% OUTPUT:
% ------- 
% - cr: the nx 1 bag of paths criticality vector containing node
%   criticalities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if A is a square matrix 
[n,m] = size(A);
if n ~= m
    display('Error : matrix A must be square');
    return;
end

%% Algorithm

% Initialization
cr = zeros(n,1); % initialize criticality vector
e  = ones(n,1); % column vector full of 1s
I  = eye(n,n); % identity matrix (same size as A)
cr = zeros(n,1); % initialize criticality vector
eps = 10e-50; % precision

d = A * e; % the outdegree matrix
Dinv = diag(1./d); % inverse of diagonal degree matrix 
    
Pref = Dinv * A; % the reference transition probabilities matrix
W = Pref .* exp(-theta*C); % % the auxiliary matrix W

Z = (I-W)\I; % the fundamental matrix

for j = 1:n % loop on nodes; compute criticality on each node in turn
    ej = I(:,j);

    zrj = Z' * ej; % row j of fundamental matrix Z viewed as a column vector
    zcj = Z * ej; % column j of fundamental matrix Z
    Z2  = Z;
    Z2(j,:) = 0; % set row j of Z to 0
    Z2(:,j) = 0; % set column j of Z to 0
    % Z2 = Z - ej*zrj' - zcj*ej' + Z(j,j)*(ej*ej'); % alternative: set row j and column j of Z to 0
    PI = Z2/sum(sum(Z2)); % the bag-of-paths probability matrix with support reduced to V \ J
    
    Zmj = Z - (zcj*zrj')/Z(j,j); % update of matrix Z when removing node j from matrix W
    PImj = Zmj/sum(sum(Zmj)); % the bag-of-paths probability matrix after removing node j
    p = PI(:); % stack columns of PI into a column vector
    pmj = PImj(:); % stack columns of PInj into a column vector

     % The criticality measure, the KL divergence between the two
     % distributions
    cr(j) = real( pmj' * log((pmj./(p + eps)) + eps) );
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
