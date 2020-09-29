function [M11, M12, M22] = getM(I,sigma)
% #### Structure Tensor ####
% To be able to detect corners in images and compute optical flow, we
% compute the structure tensor M. This function is mvg_05_01. 

% Compute (spatial) image gradients Ix and Iy using central differences.
Ix = 0.5*(I(:,[2:end end]) - I(:,[1 1:end-1]));
Iy = 0.5*(I([2:end end],:) - I([1 1:end-1],:));

% Gaussian Kernel.
k = ceil(4*sigma +1); % Kernel size given.
G = fspecial('gaussian', k, sigma);

% M is symmetric, hence M12 == M21.
M11 = conv2(Ix .* Ix, G, 'same');
M12 = conv2(Ix .* Iy, G, 'same');
M22 = conv2(Iy .* Iy, G, 'same');

end