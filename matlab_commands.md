# Matlab commands from the lecture

## Lecture 1
````matlab
%% Lecture 1
C = kron(A,B) % Kronecker product of A (m,n) and B (k,l) is (mk, nl)
Z = null(A) % Null space or kernel of matrix A
d = rank(A)
[V, D] = eig(A); % V are eigenvectors, D are eigenvalues in diag matrix
[U, S, V] = svd(A) % Singular value decomposition
X = pinv(A) % Pseudo inverse

%% Lecture 10
D = bwdist(BW) %  Assigns each pixel a number, that is the distance between that pixel and the nearest nonzero pixel of BW
```
