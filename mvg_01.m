%% #### Setup ####
% Set directory of current file as working directory
location = matlab.desktop.editor.getActiveFilename;
[path,file,ext]=fileparts(location);
cd(path);

%% #### Basic Image Processing ####
img = imread('lena.png');
disp(size(img));
imshow(img)

max = max(img, [], 'all');
min = min(img, [], 'all');

% Convert to grayscale
img_g = rgb2gray(img);

% Apply gaussian smoothing filter
rand_filter = fspecial('gaussian'); % Not just a random matrix
img_double = im2double(img_g); % Needed for convolution
img_filter = conv2(img_double, rand_filter, 'same');
imshow(img_filter);
imwrite(img_filter, 'results/img_filter.png')

% Plot with subplot
subplot(131), imshow(img),      title('Original Lena')
subplot(132), imshow(img_g), title('Grayscale Lena')
subplot(133), imshow(img_filter), title('Smoothed Lena')

%% #### Basic Operations ####
% Solve Ax = b for x
A = [2 2 0; 0 8 3];
b = [5; 15];
x = A\b; % Note the direction of the slash

B = A;
A(1,2) = 4; % Note indices start with 1, round brackets

c = 0;
for i=-4:4:4 % Start, increment, inclusive stop
  c = c + i*A.' * b; % A.' is transpose(A), both work the same
end
disp(c)

% Comparison
A .* B; % Elementwise multiplication
A' * B; % Transpose A, matrix multiplication with B

%% #### Functions - Non vectorized ####
addprimes(1,10)
approxequal([1,2;3,4], [1,2;3.1,3.5], 0.5)
addprimes_vectorized(1,10)

function sum = addprimes(s,e)
  sum = 0;
  for i=s:1:e
    if isprime(i);
      sum=sum+i;
    end
  end
end

%% #### Functions - Vectorized ####
function [out] = approxequal(x, y, eps)
    out = all(abs(x-y) < eps)
end

function [result] = addprimes_vectorized(s, e)
    z = s:e;
    result = sum(z(isprime(z)));
end
%% #### Comparison to solution ####
% Naming convention seems to be I for image.
% Don't forget to vectorize!
% Think about everything in terms of matrix multiplication.