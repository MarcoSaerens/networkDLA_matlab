function D = Alg_02_03_CommuteTimeAndEuclideanCommuteTimeDistances(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck, revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis".
%         Cambridge University Press.
%
% Description: Computes the (corrected and uncorrected) commute-time and
%              Euclidean commute-time distances between nodes of an
%              undirected and weigthed graph.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of an undirected,
%       connected, graph G.
%
% OUTPUT:
% -------
% - D : a structure containing :
%           D.CT   - a (n x n) matrix containing the average
%                    commute-time distances.
%           D.ECT  - a (n x n) matrix containing the Euclidean
%                    commute-time distances.
%           D.CCT  - a (n x n) matrix containing the corrected
%                    commute-time distances.
%           D.CECT - a (n x n) matrix containing the corrected Euclidean
%                    commute-time distances.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments

% Check if squared matrix
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared.')
end

% Check if symmetric matrix / graph is undirected
if ~isequal(A, A')
    error('The adjacency matrix is not symmetric.')
end


%% Algorithm

% Utilities
e = ones(n,1);
E = ones(n,n);
myMax = realmax;

% Diagonal matrices of degree and inverse degree
d = A*e;
Diag_d = diag(d);
Diag_d_inv = diag(1 ./ d);

% Volume of the graph
vol = sum(d);

% Laplacian matrix
L = Diag_d - A;

% Pseudoinverse of the Laplacian matrix
L_plus = ((L - (E/n))^(-1)) + (E/n);

% Average commute-time distance
diag_L_plus = diag(L_plus);
CT = vol*(diag_L_plus*e' + e*diag_L_plus' - 2*L_plus);
CT(isnan(CT)) = 0;
CT(isinf(CT)) = myMax;
D.CT = CT;

% Euclidean commute-time distance
D.ECT = CT.^(0.5);

% Corrected commute-time distance
A_norm = Diag_d_inv * A * Diag_d_inv;
diag_A_norm = diag(A_norm);
CCT = CT - vol*(Diag_d_inv*E + E*Diag_d_inv + diag_A_norm*e' + e*diag_A_norm' - 2*A_norm);
CCT(1:n+1:end) = 0;
D.CCT = CCT;

% Corrected Euclidean commute-time distance
D.CECT = CCT.^(0.5);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
