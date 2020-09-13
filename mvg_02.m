%% #### Topic ####
% Moore-Penrose pseudo-inverse
% Goal: solving the linear system Dx = b with D in R^(mx4), b in R^m,
% a vector whose components are all equal to 1, and x* = [4,-3,2,-1]^T in 
% R^4 should be one possible solution for the linear system, i.e. for any 
% row [d1, d2, d3, d4] of D: 1 = 4*d1 - 3*d2 + 2*d3 - d4
% Set of all possible solutions is given by S = {x* + v|v in kernel(D)}

%% #### Creating Data ####

% Generating a matrix D using random variables with m=4 rows
% d4 is conditioned on the others
m = 3 % number of rows
d1 = rand(m,1); % Random between 0 and 1
d2 = rand(m,1);
d3 = rand(m,1);
d4 = 4*d1 - 3*d2 + 2*d3 -1;

% Introducing small additive errors into the data
eps = 1.e-4;
d1 = d1 .* (1 + eps*rand(m,1));
d2 = d2 .* (1 + eps*rand(m,1));
d3 = d3 .* (1 + eps*rand(m,1));
d4 = d4 .* (1 + eps*rand(m,1));

D = [d1 d2 d3 d4];

%% #### Finding the coefficients x by solving the system Dx = b

% Computing the svd (singular value decomposition)
[U, S,  V] = svd(D);

% Computing the Moore-Penrose pseudo-inverse using the results of the svd
% S_dagger is Sigma_dagger from the script, it contains the inverse of the
% diagonal matrix of non-zero singular values padded with zeros to get the
% dimensions to n,m (with D being of dimensions m,n)
S_dagger = S'; % Tranpose
% .^ is elementwise potentiation
% ~= 0 is logical negation, returns 1 if element is 0, 1 otherwise
% We use the brackets with the logical negaion on both sides just to change
% the non zero values of the matrix. Would we not do this, we would get inf
% in all previously zero positions
S_dagger(S_dagger ~= 0) = S_dagger(S_dagger ~= 0) .^ -1; 
b = ones(m,1);

% Solving the least sqaures problem using pseudo inverse
D_dagger = V * S_dagger * U';
D_dagger_matlab = pinv(D);
disp(abs(D_dagger - D_dagger_matlab)) % Compare the two
% For higher m there will be lower precision


x = D_dagger * b;
disp(x)

%% #### Plotting lambdas ####
% Assume that m=3, hence we have infinately many solutions
% rank(D) = 3, dim(kernel(D)) = 1

% We know that x_min = D_daggerb is among all minimizers of ||Dx-b||^2 the
% one with the smallest norm |x|. 
% For lambda in R, x_lambda = x + lambda x one possible solution and
% e_lambda = ||Dx_lambda - b||^2 the associated error. Show that the above
% statement holds by plotting ||x_lambda|| and e_lambda.

if (m==3)
  
  v = null(D); % Null space of matrix, returns a vector v in kernel(D)
  norm(v);
  lambda = -100:1:100;
  nb_values = length(lambda);
  values_norm = zeros(nb_values,1);
  values_error = zeros(nb_values,1);
  
  for i = 1:nb_values
    x_lambda = x + lambda(i)*v;
    values_norm(i) = norm(x_lambda);
    values_error(i) = norm(D*x_lambda -b)^2;
  end
  
  % Display the graph
  figure(1)
  plot(lambda, values_norm, 'b-', lambda, values_error, 'r-');
  % Notice how all the x_lambda have the same error, and for lambda = 0, we
  % have the smallest norm possible, which we tried to show.
end

%% #### Notes ####
% To solve a linear system you would typically use x= D\b;
