%% #### Topic ####
% Epipolar lines

%% #### Data Loading and Displaying ####
% Read images:
I0 = double(imread('batinria0.pgm'));
I1 = double(imread('batinria1.pgm'));
[w,h] = size(I1);

% Parameters from the file 'calibration.txt'. Rotation/Translation
% transform points form camera coordinates to some global coordinate frame.
% Camera parameters left image:
K0 = [844.310547 0 243.413315; 0 1202.508301 281.529236; 0 0 1];
R0 = [0.655133, 0.031153, 0.754871;0.003613, 0.999009, -0.044364;-0.755505, 0.031792, 0.654371];
T0 = [-793.848328; 269.264465; -744.572876];
% Camera parameters right image:
K1 = [852.721008 0 252.021805; 0 1215.657349 288.587189; 0 0 1];
R1 = [0.739514, 0.034059, 0.672279;-0.006453, 0.999032, -0.043515;-0.673111, 0.027841, 0.739017];
T1 = [-631.052917; 270.192749; -935.050842];

% Show the two images and select a point in  image zero. Using ginput to
% retrieve the image coordinates of a mouse click.
% Image 0.
figure; imshow(uint8(I0));
hold on;
[x0,y0] = ginput(1);
plot(x0,y0,'r+');
hold off;

% Image 1.
% figure; imshow(uint8(I1));

%% #### Compute the Epipolar Line l_1 in the Second Image ####
% Goal: Computing l = F*x_0.

% Fundamental matrix F:
g0 = [R0 T0; 0 0 0 1];
g1 = [R1 T1; 0 0 0 1];
g = inv(g1) * g0;
T = g(1:3,4);
R = g(1:3,1:3);
F = inv(K1)' * hat(T) * R * inv(K0);

% Epipolar line for x_0:
l = F * [x0; y0; 1];

% Draw epipolar line on I1:
figure; imshow(uint8(I1));
hold on;
m = -l(1)/l(2);
b = -l(3)/l(2);
y0 = m * 1 + b;
y1 = m * w + b;
line([1 w],[y0 y1]);

% TODO: solve bonus question.

%% #### Functions ####
% Create a skew-symmetric matrix.
function A = hat(v)
A = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0];
end