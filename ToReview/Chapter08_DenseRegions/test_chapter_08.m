%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Masashi Shimbo and Marco Saerens (2018).
% 
% 
% Description: Tests for all code in chapter 04.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of elements for tests
n = 10;

% utilities
e = ones(n, 1); % vector of ones
I = eye(n,n); % identity matrix
eps = 10^(-100); % precision

rng(72) % set seed of random generator

% % An adjacency matrix for a connected unweighted, undirected, graph with n nodes
% A = rand(n,n);
% A = A > 0.6325; % sparsify & binarize
% A(1:n+1:end) = 0; % nullify main diagonal (remove self-loops)
% A(2:n+1:end) = 1; % ensure connectedness
% A = double(A | A'); % symmetrize A

% % An adjacency matrix for a connected weighted undirected graph with n nodes
% A = rand(n,n);
% A = (A > 0.6325) .* A; % sparsify
% A(1:n+1:end) = 0; % nullify main diagonal (remove self-loops)
% A(2:n+1:end) = 1; % ensure connectedness
% A = (A + A')/2; % symmetrize A

% Some other test examples

% A = [ 0 1 1 0 0 0
%       1 0 1 0 0 0
%       1 1 0 1 0 0
%       0 0 1 0 1 1
%       0 0 0 1 0 1
%       0 0 0 1 1 0 ];

% Example of the book for k-cores
A = [ 0 1 0 0 0 0 0 0 0 0 0 0   % 01
      1 0 1 0 0 1 0 0 0 0 0 0   % 02
      0 1 0 0 0 0 0 0 0 0 0 0   % 03
      0 0 0 0 1 0 0 0 0 1 0 0   % 04
      0 0 0 1 0 1 1 1 0 0 0 0   % 05
      0 1 0 0 1 0 1 1 0 0 0 0   % 06
      0 0 0 0 1 1 0 1 0 1 1 0   % 07
      0 0 0 0 1 1 1 0 1 0 1 0   % 08
      0 0 0 0 0 0 0 1 0 0 0 0   % 09
      0 0 0 1 0 0 1 0 0 0 0 0   % 10
      0 0 0 0 0 0 1 1 0 0 0 1   % 11
      0 0 0 0 0 0 0 0 0 0 1 0 ];% 12

%% Test 3 Alg_08_03_KCore(A)

k = 3;
u = Alg_08_03_KCore(A, k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
