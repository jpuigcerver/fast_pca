function [Z, Zh, V, D] = pca(X)
% Mean and standard deviation of the data
m = mean(X);
s = std(X);
% Zero-Mean centering
B = bsxfun(@minus, X, m);
% Normalization of the mean-centered data
Bh = bsxfun(@rdivide, B, s);
% Covariance matrix of the zero-mean data
C = cov(B);
% Eigenvectors and eigenvalues of the cov. matrix
[V,D] = eig(C);
D = diag(D);
% Sort eigenvalues and eigenvectors
D = D(end:-1:1);
V = fliplr(V);
% Project data
Z = B * V;
% Project normalized data
Zh = Bh * V;
end
