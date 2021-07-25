# Matlab commands from the lecture

## Lectures
```matlab
%% Lecture 1
C = kron(A,B) % Kronecker product of A (m,n) and B (k,l) is (mk, nl)
Z = null(A) % Null space or kernel of matrix A
d = rank(A)
[V, D] = eig(A); % V are eigenvectors, D are eigenvalues in diag matrix
[U, S, V] = svd(A) % Singular value decomposition, sizes for A: (m,n), U(m,m), S(m,n), V(n,n)
X = pinv(A) % Pseudo inverse

%% Lecture 10
D = bwdist(BW) %  Assigns each pixel a number, that is the distance between that pixel and the nearest nonzero pixel of BW
```

## Tutorials
```matlab
x = A\b % Solve system of linear equations Ax = b for x
x = pinv(A)*b % Same way to solve
x = inv(A'*A) * A' * b % Explicit way to construct the pseudo inverse for a (n,m) matrix

T = [eye(3) [-0.5 -0.5 1]'];
```

## Exams
```matlab
y = A\b % [U, S, V] = svd(A); x = V*pinv(S)*U'*b is equal up to numerical percision, A is (3,2), b is (3,1)
```
## Other operands
```matlab
B = A .* A % Elementwise multiplication, ./ for elementwise division, .^ for potentiation
A_trans = A' % Transpose
A_inv = inv(A)
A_stack = A(:) % Stacks columns
T = g(1:3,4) % Indexing, extracting T out of g (row 1-3, col 4).
Z = cross(X,Y) % Crossproduct, creates new vector if X, Y are 3x1 vectors
scalar_product = dot(A,B) % Equivalent to A'*B

zeros(n) % (n,n) matrix of zeros
zeros(m,n) % (m,n) matrix of zeros
zeros(size(A)) % Size like a but all zeros
expm(X) % Matrix exponential of X, Although it is not computed this way, if X has a full set of eigenvectors V with corresponding eigenvalues D, then [V,D] = eig(X) and expm(X) = V*diag(exp(diag(D)))/V. expm(X) yields triangular matrix, diagonal elements are same as for exp(X), other entries are different


trace(A) % Trace
B_stacked = B(:) % Stacks the columns of B
rank(A)
inv(A)
eye(n,m) % Identity matrix of this size
ones(n,m)
magic(n) % N-by-N matrix constructed from the integers 1 through N^2 with equal row, column, and diagonal sums. Produces valid magic squares for all N > 0 except N = 2.
sum(A) % This returns column sums, e.g. A=[1,2;3,4]; sum(A) -> [4, 6]
sum(sum(A)) % For square matrix this is sum of all values
% Multiline matrix initialization
K = [540 0 320; ...
     0 540 240; ...
     0   0   1];

w_hat = hat(w) % Use the function from below.
h = @ hat; % Create a function handle, can do h(w) instead of hat(w) now.
g_ATAN_1 = @(r) (1./(w1*r) .* atan(2*tan(w1/2)*r)); % Inline function.
result = g_ATAN_1(40) % This is how to call the inline function.



% Create skew-symmetric matrix. 
function A = hat(v) 
A = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0];
end

% Special indexing, useful for gradients.
A = [0 1 2 3 4 5 6 7 8 9; 10 11 12 13 14 15 16 17 18 19]
A(1,[2:end]) == [1 2 3 4 5 6 7 8 9] % Starting at 2nd until last.
A(1,[2:end end]) == [1 2 3 4 5 6 7 8 9 9] % Repeating last element, used for gradient calculation to enlarge image.

% Calculate image gradients
Ix = 0.5*(I(:,[2:end end]) - I(:,[1 1:end-1]));
Iy = 0.5*(I([2:end end],:) - I([1 1:end-1],:));

% .^ is elementwise potentiation
% ~= 0 is logical negation, returns 1 if element is 0, 1 otherwise.
% We use the brackets with the logical negaion on both sides just to change
% the non zero values of the matrix. Would we not do this, we would get inf
% in all previously zero positions
S_dagger(S_dagger ~= 0) = S_dagger(S_dagger ~= 0) .^ -1;
```

## General Remarks

### Vectorization
Avoid for loops and make use of matlabs vector based operations.

### Order of operations
Think about the order of operations to decrease the amount of operations needed to solve a problem. 

Example: let X be of size (100000,3), w_hat be a (3,3) skew-symmetric matrix and R be a (3,3) rotation matrix.
`X * (w_hat * R)'` is about twice as fast as `(X * R') * w_hat'` just due to the order of operators. The first example only does once a multiplication with a large matrix of size (100000,3), compared to the second one, that does it twice, since its left bracket is a large matrix again.