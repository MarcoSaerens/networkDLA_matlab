%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Author: Guillaume Guex (2017).
% 
% 
% Description: Tests for all code in chapter 7.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of elements for tests
n = 100;

% utilities
e = ones(n, 1);
I = eye(n);
H = I - e*e'/n;

%% Test 1

% Distance matrix
X = rand(n, 2);
D = zeros(n);
for i = 1:n
    for j = 1:n
        D(i,j) = (X(i, :) - X(j, :)) * (X(i, :) - X(j, :))';
    end
end
% number of groups
m = 10;

% TEST Alg_07_01_StandardkMeansClustering(D, m)
groups = Alg_07_01_StandardkMeansClustering(D, m);

% Display
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
scatter(X(:,1), X(:,2),[], node_color)

%% Test 2

K = - 1/2 * H * D * H;
K = 1/2 * (K + K');

% TEST Alg_07_02_KernelkMeansClustering(K, m)
groups = Alg_07_02_KernelkMeansClustering(K, m);

% Display
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 3

% TEST Alg_07_03_IterativeKernelkMeansClustering(K, m)
groups = Alg_07_03_IterativeKernelkMeansClustering(K, m);

% Display
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 4

% Neighbours graph
n_neighbours = 5;
A = zeros(n, n);
for i = 1:n
   [~, ind] = sort(D(i, :));
   A(i, :) = (ind < n_neighbours + 2);
end
A = (A + A');
A(1:(n+1):end) = 0;
A = (A > 0) ./ (D + 10);

% lambda
lambda = 0;

% TEST Alg_07_04_LabelPropagationClustering(A, lambda)
groups = Alg_07_04_LabelPropagationClustering(A, lambda);

% Display
[~, m] = size(groups);
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 5

% TEST Alg_07_05_MarkovClustering(A, q, r)
groups = Alg_07_05_MarkovClustering(A, 2, 2);

% Display
[~, m] = size(groups);
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 6

% TEST Alg_07_06_RegularizedMarkovClustering(A, r)
groups = Alg_07_06_RegularizedMarkovClustering(A, 1.5);

% Display
[~, m] = size(groups);
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 7

% Creates U_init
u_1 = (X(:,1) > 0.5);
u_2 = 1 - u_1;
U_init = [u_1, u_2];

% TEST Alg_07_07_KernighanLinClustering(A, U_init)
groups = Alg_07_07_KernighanLinClustering(A, U_init);

% Display
[~, m] = size(groups);
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 8

% TEST Alg_07_08_RatioCutSpectralClustering(A, m, option)
groups = Alg_07_08_RatioCutSpectralClustering(A, 5);

% Display
[~, m] = size(groups);
group_color = linspace(1,10,m);
node_color = zeros(1, n);
for g = 1:m
    node_color(logical(groups(:, g))) = group_color(g);
end
figure
scatter(X(:,1), X(:,2),[], node_color)

%% Test 9

%% Test 10

%% Test 11

% TEST 
Alg_07_11_LatentClassClustering(A, m)
