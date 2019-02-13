function st = Alg_02_10_MarkovDiffusionSquareDistance(A, w, t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Maxime Duyck revised by Guillaume Guex (2017).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis".
%         Cambridge University Press.
%
% Description: Computes the Markov diffusion square distance and kernel
%              matrix of a weighted directed graph.
%
% INPUT:
% -------
% - A : the (n x n) weighted adjacency matrix of a directed graph G.
% - w : the weighting factors w_k >= 0 associated with each node k in
%       vector w. When the Markov chain is regular, the usual weighting
%       factor is w_k = 1/pi_k where pi_k is entry k of the stationary
%       distribution.
% - t : the considered number of steps (or transitions) of the Markov
%       process.
%
% OUTPUT:
% -------
% - st : a structure containing :
%           st.D_MD  - the (n x n) Markov diffusion squared distances
%                      matrix, original version.
%           st.D_MDA - the (n x n) Markov diffusion squared distances
%                      matrix, time-averaged version.
%           st.K_MD  - the (n x n) Markov diffusion kernel matrix,
%                      original version.
%           st.K_MDA - the (n x n) Markov diffusion kernel matrix,
%                      time-averaged version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks of arguments

% Check if squared matrix
[n, m] = size(A);
if n ~= m
    error('The adjacency matrix is not squared')
end

% Check if length of w is correct
if length(w) ~= n
    error('The weigthing vector has not the correct dimension')
end

% Check if the weighting factors have correct values
if min(w) < 0
    error('The weighting vector has at least a negative value')
end

% Check if t is >= 1
if t < 1
    error('The number of step should be >= 1')
end

%% Algorithm

% Utilities
e = ones(n, 1);
Diag_d = diag(A*e);
I = eye(n);

% The transition matrix associated with A
P = Diag_d^(-1) * A;

% The centering matrix
H = I - (e*e')/n;

% The t-order Neumann series of A, averaged by t
Z_t = zeros(n);
P_pow_t = eye(n);
for i = 1:t
    P_pow_t = P*P_pow_t;
    Z_t = Z_t + P_pow_t;
end
Z_t = Z_t / t;

% The diagonal matrix containing the weights
Diag_w = diag(w);

% The original version of Markov diffusion squared distances matrix
Pt_Dw_Pt = P_pow_t * Diag_w * P_pow_t';
st.D_MD = diag(Pt_Dw_Pt)*e' + e*diag(Pt_Dw_Pt)' - 2*Pt_Dw_Pt;

% The time-averaged version of Markov diffusion squared distances matrix
Zt_Dw_Zt = Z_t * Diag_w * Z_t';
st.D_MDA = diag(Zt_Dw_Zt)*e' + e*diag(Zt_Dw_Zt)' - 2*Zt_Dw_Zt;

% The original version of centered Markov diffusion kernel matrix
st.K_MD = H * Pt_Dw_Pt * H;

% The time-averaged version of centered Markov diffusion kernel matrix
st.K_MDA = H * Zt_Dw_Zt * H;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%