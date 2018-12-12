function [Zx, Zy] = Alg_09_01_CorrespondenceAnalysis(W, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Francois Fouss, revised by Marco Saerens (2018).
%
% Source: Francois Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes a simple correspondence analysis of the 
% bi-adjacency matrix W of an undirected, weighted, bipartite graph.
% the output is an embedding of the left and right nodes in a p-dimensional
% space
%
% INPUT:
% ------ 
% - W: the (nx x ny) weighted biadjacency matrix of an undirected bipatite
%   graph G, assumed strongly connected. x denotes left nodes while y
%   denotes right nodes.
% - p: the dimensionality of the embedding space.
% 
% OUTPUT:
% -------
% - Zx, Zy: matrices (nx x p and ny x p) containing the coordinates of the
%   left set of nodes and the right set of nodes in the embedding space.
%   Each line is the coordinate vector of the corresponding node. The
%   columns correspond to the different dimensions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Algorithm
[nx,ny] = size(W);
 
ex = ones(nx,1); % column vector full of 1s of the appropriate size
ey = ones(ny,1); % column vector full of 1s of the appropriate size
Zx = []; % will contain the coordinates of left nodes
Zy = []; % will contain the coordinates of right nodes
 
N = ex' * W * ey; % volume of interactions
dx = W  * ey; % outdegree vector of left (x) nodes
dy = W' * ex; % outdegree vector of right (y) nodes

Dx = diag(dx); % diagonal matrix with outdegrees of left nodes
Dy = diag(dy); % diagonal matrix with outdegrees of right nodes
 
Dx_inv = diag(1./dx); % diagonal matrix with inverse outdegrees of left nodes
Dy_inv = diag(1./dy); % diagonal matrix with inverse outdegrees of right nodes
 
Px = Dx_inv * W  * Dy_inv * W'; % the transition matrix between left nodes
Py = Dy_inv * W' * Dx_inv * W; % the transition matrix between right nodes

% Compute p right dominant eigenvectors of Px and sort them in decreasing order
[Vx,Dvx] = eigs(Px,p+1);
[value,ind] = sort(diag(Dvx),'descend');
Dux = diag(value);
Ux  = Vx(:,ind); % Ux contains the eigenvectors as columns

% Compute p right dominant eigenvectors of Py and sort them in decreasing order
[Vy,Dvy] = eigs(Py,p+1);
[value,ind] = sort(diag(Dvy),'descend');
Duy = diag(value);
Uy  = Vy(:,ind); % Uy contains the eigenvectors as columns

for k = 2:(p+1)
  ux_k = Ux(:,k);
  uy_k = Uy(:,k);
  
  % normalize the eigenvectors
  ux_k = sqrt(Dux(k,k)) * ux_k/sqrt((ux_k' * Dx * ux_k)/N);
  uy_k = sqrt(Duy(k,k)) * uy_k/sqrt((uy_k' * Dy * uy_k)/N);
  
  % stack the p eigenvectors in the Z matrices
  Zx = [Zx,ux_k];
  Zy = [Zy,uy_k];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
