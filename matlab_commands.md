# Matlab commands from the lecture

## Lectures
```matlab
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

## Tutorials
```matlab
x = A\b % Solve system of linear equations Ax = b for x
x = pinv(A)*b % Same way to solve
x = inv(A'*A) * A' * b % Explicit way to construct the pseudo inverse for a (n,m) matrix

T = [eye(3) [-0.5 -0.5 1]'];
```

## Other operands
```matlab
B = A .* A % Elementwise multiplication, ./ for elementwise division
A_trans = A' % Transpose
A_inv = inv(A)
A_stack = A(:) % Stacks columns
T = g(1:3,4) % Indexing, extracting T out of g (row 1-3, col 4).

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
```

## General Remarks

### Vectorization
Avoid for loops and make use of matlabs vector based operations.

### Order of operations
Think about the order of operations to decrease the amount of operations needed to solve a problem. 

Example: let X be of size (100000,3), w_hat be a (3,3) skew-symmetric matrix and R be a (3,3) rotation matrix.
`X * (w_hat * R)'` is about twice as fast as `(X * R') * w_hat'` just due to the order of operators. The first example only does once a multiplication with a large matrix of size (100000,3), compared to the second one, that does it twice, since its left bracket is a large matrix again.