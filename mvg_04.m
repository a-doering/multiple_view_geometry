%% #### Topic ####
% This exercise covers image formation, radial distortion and image
% rectification.

%% #### Image Formation ####

% a) ----------------------------
% Load data using the provided openOFF function. V are vertices, F are
% faces.
[V, F] = openOFF('model.off');
% Number of vertices:
N_vert = size(V,1);

% Translate the model vertices by vector T using homogeneous coordinates
% and a 3x4 rigid body transformation
T = [eye(3) [-0.5 -0.5 1]'];
% Transpose for easier vector-wise multiplication
V_hom = [V'; ones(1, N_vert)]; % Size (4, 19105)
V_shifted = T*V_hom; % Shift and perspective projection, cut 4th coordinate

% b) ----------------------------
% Computing the perspective projection of the new model by a camera with
% intrinsic parameters and visualized projection.

% Camera intrinsic parameters:
K = [540 0 320; ...
     0 540 240; ...
     0   0   1];

% Apply camera matrix and convert from homogeneous to pixel coordinates.
V_proj_hom = K*V_shifted;
% Take the x,y values and divide by the respective z values
V_proj_uv = V_proj_hom(1:2,:) ./ V_proj_hom(3,:);

% Visualization.
figure(1)
subplot(121)
patch('Vertices', V_proj_uv', 'Faces', F)
% To get an image that corresponds to our convention of the center of the 
% top-left pixel having coordinate (0,0), we set the viewport (axis) to a 
% range that starts at -0.5.
axis equal, axis([0 640 0 480]-0.5) % Dimensions were given.
title('Perspective projection')

% c) ----------------------------
% Parallel projection model, the generic perspective projection pi is
% replaced by a parallel projection orthogonal to the z-axis. Same as in
% b), just other projection.

% Do parallel projection by setting z to 1, applying camera intrinsics and
% converting to pixel coordinates.
V_proj_parallel = [V_shifted(1:2,:); ones(1,N_vert)];
V_proj_parallel_hom = K * V_proj_parallel;
V_proj_parallel_uv = V_proj_parallel_hom(1:2,:) ./ V_proj_parallel_hom(3,:);

% Visualization.
subplot(122)
patch('Vertices', V_proj_parallel_uv', 'Faces', F)
axis equal, axis([0 640 0 480]-0.5)
title('Parallel projection')

%% #### Radial Distortion and Image Rectification ####

% We compute a rectified image from a radially distorted image. Given a
% distorted image Id with projection function pi_d as defined before and
% camera intrinsics Kd this amounts to computing a virtual, undistorted
% image Inew with generic projection function pi and a given camera matrix
% K_new.

% a) ----------------------------

% Load image and provided camera intrinsics and visualize the image.
Id1 = imreadbw('img1.jpg');
Kd1 = [388.6 0 343.7; ...
       0 389.4 234.6; ...
       0     0     1];
w1 = 0.92646; % Omega from the ATAN model (g = g_ATAN) from the theory part.

figure(2)
subplot(121)
% NOTE: imagesc(Id1) in MATLAB will plot the image in the range 
%       (0.5, width+0.5) x (0.5, height+0.5) such that the center of the
%       top-left pixel is at (1,1). For visualization this is not a
%       problem, but later, when we lookup/interpolate pixel values, we
%       need to account for this convetion vs (0,0) for the top-left pixel.
imagesc(Id1), axis image, colormap gray
title('Distorted image')
% Marvel at how the distortion curves straight lines.

% b) ----------------------------
% Compute a virtual rectified image Inew of dim 1024x768 with a projection
% function according to a pinhole camera model and intrinsic parameters
% K_new.

% Define the distortion function inline (use 'res = g_ATAN_1(val)') with
% parameter w1. Can be applied to a whole vector of r-values.
g_ATAN_1 = @(r) (1./(w1*r) .* atan(2*tan(w1/2)*r));

% Desired camera instrinsics for rectified image (given).
K_new = [250 0 512; ...
         0 250 384; ...
         0   0   1];
       
% Create meshgrid starting with 0,0 in the top left to have a place for all
% pixels in the image we want to create. For each we need to later lookup
% an (interpolated) intensity value in the distorted image. Create vector
% uv_hom to have all pixels in homogeneous coordinates in one vector for
% easier manipulation.

% u is a matrix (768x1024) where each row is a copy of 0:1023, v also a 
% matrix, but each column is a copy of 0:767.
[u,v] = meshgrid(0:1023, 0:767); 
N_img = 1024*768;
uv_hom = [u(:) v(:) ones(N_img, 1)]';


% Unproject image coordinates of ideal pinhole camera to generic image
% plane (at Z=1). 
X_generic = inv(K_new) * uv_hom;

% Compute norm of the undistorted image coordinates.
r = sqrt(X_generic(1,:).^2 + X_generic(2,:).^2);

% Apply distortion, ignore z coordinate of X_generic, since it is 1 for
% all.
X_d1 = [g_ATAN_1(r) .* X_generic(1:2,:); ones(1, N_img)];

% Project distorted coordinates to actual image.
uv_d1_hom = Kd1 * X_d1;

% Find pixel values for each point in uv_d1_hom by linear interpolation.
% Need to ensure that top-left corner has coordinates (0,0). Also we ignore
% the z cordinate of the homogeneous vectors in uv_d1_hom, since we know
% they will be 1. Pixels outside of original image will be set to black
% (0). Reshape vector of pixel values to rectangualr image.
[Hd1, Wd1] = size(Id1);
[grid_u_d1, grid_v_d1] = meshgrid(0:Wd1-1, 0:Hd1-1);
Inew = interp2(grid_u_d1, grid_v_d1, Id1, uv_d1_hom(1,:), uv_d1_hom(2,:), 'linear', 0);
Inew = reshape(Inew, size(u));

% Visualizing the rectified image and save it to disc.
subplot(122)
imagesc(Inew), axis image, colormap gray
title('Undistorted image')
imwrite(Inew,'results/img1_undist.jpg')


% c) ----------------------------
% Same procedure as b, but the distortion is modeled by a polynomial, we
% have different camera intrinsics and work on image2.

Id2 = imreadbw('img2.jpg');
Kd2 = [279.7 0 347.3; ...
       0 279.7 235.0; ...
       0     0     1];

% Polinomial distortion function, working on a vector of r values.
g_pol_2 = @(r) 1 - 0.3407*r + 0.057*r.^2 - 0.0046*r.^3 + 0.00014*r.^4;
   
figure(3)
subplot(121)
imagesc(Id2), axis image, colormap gray
title('Distorted image')

% Apply distortion, ignore z coordinate of X_generic, since it is 1 for
% all.
X_d2 = [g_pol_2(r) .* X_generic(1:2,:); ones(1, N_img)];

% Project the distorted coordinates to the actual image.
uv_d2_hom = Kd2 * X_d2;

% Find pixel values for each point in uv_d2_hom by linear interpolation.
% Need to ensure that top-left corner has coordinates (0,0). Also we ignore
% the z cordinate of the homogeneous vectors in uv_d2_hom, since we know
% they will be 1. Pixels outside of original image will be set to black
% (0). Reshape vector of pixel values to rectangualr image.
[Hd2, Wd2] = size(Id2);
[grid_u_d2, grid_v_d2] = meshgrid(0:Wd2-1, 0:Hd2-1);
Inew2 = interp2(grid_u_d2, grid_v_d2, Id2, uv_d2_hom(1,:), uv_d2_hom(2,:), 'linear', 0);
Inew2 = reshape(Inew2, size(u));

% Visualize the rectified image and save to disc.
subplot(122)
imagesc(Inew2), axis image, colormap gray
title('Undistorted image')
imwrite(Inew2,'results/img2_undist.jpg')

% TODO: bonus