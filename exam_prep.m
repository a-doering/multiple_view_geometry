%% Epipolar constraint, if res_1 is all zero --> then x1, x2 are images of each other
T = [5; 4; 1];
R = [30;20;1];
% Replace values in x_2, x_1
x_2 = [20; 68/5; 1];
x_1 = [-16; 25; -20];

res_1 = x_2' * cross(T, R) * x_1
%%
A = [  ]
B = [ ]

%%
g = [0 0 1 4; -1 0 0 2; 0 -1 0 3; 0 0 0 1]

%% 
m = 5
n = 7
P = kron(eye(n), eye(m))
A = rand(m,n)
res = kron(A(:), A(:))'*P(:)

%res == sum(sum(A*A'))%trace(A*A')
trace(A*A')
trace(A*A')
%%
[x,y] = eig(rand(3))

%% New Exam Prep

A = [1,2;3,4]
diag(A)

format rational
%% testing
a=[0.65 0.75]
