%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Masashi Shimbo and Marco Saerens (2018).
% 
% 
% Description: Tests for all code in chapter 09.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% W = [0 1 1 0 0 0 0;
%      1 0 1 0 0 0 0;
%      1 1 0 1 0 0 0;
%      0 0 1 0 1 0 0;
%      0 0 0 1 0 1 1;
%      0 0 0 0 1 0 1;
%      0 0 0 0 1 1 0]; % example of a symmetric matrix with 2 communities

%  W = [1 0 0 0;
%       0 1 1 0;
%       0 0 1 1;
%       1 1 1 0;
%       0 0 1 1;
%       0 0 0 1]; % example of the book, page 391
 
%  W = [0   6   1   0;
%       0   9   1   0;
%       1  12   6   0;
%       7  11   8   1;
%       8   7  23   3;
%       8   6   7   7;
%       19  1   5  14;
%       8   0   1  16;
%       1   0   0  11]; % the well-known "cheese" dataset for testing correspondence analysis

% % Number of elements
% [nx,ny] = size(W);

% Utilities
eps = 10^(-100); % precision

% A random biadjacency matrix for a connected weighted undirected graph
% with nx + ny nodes

rng(72) % set seed of random generator
nx = 10; ny = 20;

W = rand(nx, ny);
W(W <= 0.6325) = 0; % sparsify
% 
% C = 1./(W + eps); % compute cost matrix

[Zx,Zy] = Alg_09_01_CorrespondenceAnalysis(W,3)
