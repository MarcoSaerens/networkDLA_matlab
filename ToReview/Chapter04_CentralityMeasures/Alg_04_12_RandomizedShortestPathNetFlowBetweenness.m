function bet = Alg_04_12_RandomizedShortestPathNetFlowBetweenness(Pref, C, theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Ilkka Kivimaki (2015), revised by Marco Saerens (2018).
%
% Source: François Fouss, Marco Saerens and Masashi Shimbo (2016).
%         "Algorithms and models for network data and link analysis". 
%         Cambridge University Press.
%
% Description: Computes the net flow RSP betweenness containing
% the net flow through each intermediate node, accumulated over
% all source-destination pairs.
%
% INPUT:
% ------ 
% - Pref: the (n x n) transition probability matrix of a weighted,
%   strongly connected, undirected, graph G.
% - C: the (n x n) cost matrix C associated to G (default = 1./A).
% - theta: the (scalar) non-negative inverse temperature parameter.
% 
% OUTPUT:
% -------
% - bet: the RSP net flow betweenness vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some checks 
[m,n] = size(Pref);
if m ~= n
    error('Matrix Pref must be square')
end

%% Algorithm
W = Pref .* exp(-theta*C); % RSP's W
I = eye(n); % identity matrix (same size as W)
e = ones(n,1);
bet = zeros(n,1);

Z = (I - W)\I; % the fundamental matrix Z = inv(I - W)
  
if any(any(Z<=0 | isinf(Z)))
    error('Z contains zero or Inf values - either graph is not strongly connected or beta is too small/large');
end

Zdiv = 1./Z; % the matrix containing elements 1/Z(i,j)
dZdiv = diag(Zdiv); % column vector containing the diagonal elements 1/Z(i,i)
                          
% Compute betweenness for each node; we loop over all edges j->jj (j->j' in the
% book)
for jj = 1:n
    for j = find(Pref(jj,:)) % loop on neighbors j of jj (non-zero transition probabilities)
        % Flow in edge j->jj and jj->j for all source-destination nodes
        Nj2jj = W(j,jj) * ( ( (Z(:,j)*Z(jj,:)) .* Zdiv ) - ( e * (Z(:,j).*(Z(jj,:))'.*dZdiv)' ) ); % flow in edge j->jj
        Njj2j = W(jj,j) * ( ( (Z(:,jj)*Z(j,:)) .* Zdiv ) - ( e * (Z(:,jj).*(Z(j,:))'.*dZdiv)' ) ); % flow in edge jj->j
        
        N = abs(Nj2jj - Njj2j); % net flow in edge i->jj for all source-destinations
        
        bet(jj) = bet(jj) + (e - I(:,jj))' * N * (e - I(:,jj)); % update betweenness
    end
end
         
% The RSP net flow betweenness vector
bet = (1/(2*(n-1)*(n-2))) * bet;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%