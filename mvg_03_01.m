%% #### Topic 1 ####
% Homogeneous matrices
% Goal: Rotate model of vertices around its center (i.e. the mean of its 
% vertices) for rotation angles alpha, beta, gamma around the x,y,z axis. 
% Use homogeneous cordinates and describe the overall transformation by a 
% single matrix.

%% #### Load Model ####
% Function openOFF is provided in another file.
[V,F] = openOFF('model.off', '');

%% #### Rotate/Translate Examples ####
% Display model
close all;
figure
display(V, F)

% Rotate first 45° around x, then 120° around z
figure
V1 = transform(V, [deg2rad(45) 0 0]', [0 0 0 ]');
V1 = transform(V, [0 0 deg2rad(120)]', [0 0 0]');
display(V1, F)

% Now first around z, then x.
figure
V2 = transform(V,  [0           0 deg2rad(120)]', [0 0 0]');
V2 = transform(V2, [deg2rad(45) 0 0           ]', [0 0 0]');
display(V2, F)
% --> This shows that the order in which the rotations are applied matters.

% Translate.
figure
V3 = transform(V, [0 0 0]', [0.5 0.2 0.1]');
display(V3, F)

%% #### Display Model ####
function display(V,F)
    C = 0.3*ones(size(V,1),3);
    patch('Vertices', V, 'Faces', F, 'FaceVertexCData', C);
    axis equal;
    shading interp;
    camlight right;
    camlight left;
end

%% #### Transforms ####
function V = transform(V, angles, t)
  % Rotate vertices V around their center with provided angles around the
  % x,y,z axis and then translate by t

    centroid = mean(V)';
    t_zero = zeros(3, 1);

    % Right to left: shift to center, rotate by x, y, z, shift back,
    % translate by translation.
    T = SE3(eye(3), t) * ...
        SE3(eye(3), centroid) * ... 
        SE3(rotmatz(angles(3)), t_zero) * ...
        SE3(rotmaty(angles(2)), t_zero) * ...
        SE3(rotmatx(angles(1)), t_zero) * ...
        SE3(eye(3), -centroid);
    
    Vhom = [V ones(size(V, 1), 1)]';
    Vhom = T * Vhom;
    
    % We can simply drop the last element since it will still be 1 in this
    % case. In general to convert to euclidian coordinates we need to
    % divide by the last element.
    V = Vhom(1:3,:)';
end

function T = SE3(R, t)
    T = eye(4);
    T(1:3, 1:3) = R;
    T(1:3, 4) = t;
end

function R = rotmatx(a)
    R = [1      0       0;
         0 cos(a) -sin(a);
         0 sin(a)  cos(a)];
end

function R = rotmaty(a)
    R = [ cos(a) 0 sin(a);
               0 1      0;
         -sin(a) 0 cos(a)];
end

function R = rotmatz(a)
    R = [cos(a) -sin(a) 0;
         sin(a)  cos(a) 0;
              0       0 1];
end