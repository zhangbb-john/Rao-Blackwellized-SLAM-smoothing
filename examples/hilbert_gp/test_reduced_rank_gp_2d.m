clear;
close all;
% Parameters
m = 128;         % Number of eigenfunctions
d = 2;          % Dimensionality
L = [4, 4];     % Domain size [-L1, L1] x [-L2, L2]
x = [3.5, 3.5]; % Test point 1
x_prime = [0.5, 0.8]; % Test point 2

% Generate eigenvalues and eigenfunctions
[eigenval, eigenfun, ~, NN] = domain_cartesian_dx(m, d, L);

% The eigenvalues (unscaled)
lambda_raw = eigenval(NN);

% Define the spectral density function for SE kernel
sigma2 = 1; % Signal variance
lengthScale = 1;    % Length-scale
% S = @(w, lengthScale, magnSigma2) ...
%     magnSigma2 * (2*pi*lengthScale^2)^(d/2) * exp(-w.^2 * lengthScale^2 / 2);
S = @(w,lengthScale,magnSigma2) ...
    magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
% k = S(sqrt(lambda),lengthScale,magnSigma2);
% Scale eigenvalues with spectral density
k = zeros(size(lambda_raw)); % Initialize scaled eigenvalues
for n = 1:m
    w = norm(pi * NN(n, :) ./ (2 * L)); % Frequency w corresponding to index NN
    k(n) = S(sqrt(lambda_raw(n)), lengthScale, sigma2); % Scale eigenvalue by spectral density
end

% Compute the approximate SE kernel value
kernel_approx = 0; % Initialize kernel approximation
for n = 1:m
    phi_x = eigenfun(NN(n, :), x);        % Eigenfunction at x
    phi_x_prime = eigenfun(NN(n, :), x_prime); % Eigenfunction at x'
	disp(['phi_x is ', num2str(phi_x), '; phi_x_prime is ', num2str(phi_x_prime),...
		'; (k(n)) * phi_x * phi_x_prime is ', num2str((k(n)) * phi_x * phi_x_prime)]);
    kernel_approx = kernel_approx + k(n) * phi_x * phi_x_prime;
end

% Compute the true SE kernel value for comparison
dist_sq = sum((x - x_prime).^2); % Squared distance
kernel_true = sigma2 * exp(-dist_sq / (2 * lengthScale^2)); % True SE kernel

% Display results
disp(['Approximated SE kernel value: ', num2str(kernel_approx)]);
disp(['True SE kernel value: ', num2str(kernel_true)]);

%% Visualization
% Create a grid of points for visualization
n_points = 20;
x_grid = linspace(-L(1), L(1), n_points);
y_grid = linspace(-L(2), L(2), n_points);
[X, Y] = meshgrid(x_grid, y_grid);

% Compute approximate SE kernel on the grid
K_approx = zeros(n_points, n_points);
K_true = zeros(n_points, n_points);

for i = 1:n_points
    for j = 1:n_points
        x1 = [X(i, j), Y(i, j)];
		disp(['point is [', num2str(X(i, j)), ', ', num2str(Y(i, j)), ']']);
		
        % Approximate kernel
        k_approx = 0;
        for n = 1:m
            phi_x1 = eigenfun(NN(n, :), x1);% nearby points
            phi_x = eigenfun(NN(n, :), x); % Fixing the reference point x
            k_approx = k_approx + k(n) * phi_x1 * phi_x;
			disp([' This is No. ', num2str(n), ' basis function: phi_x1 * phi_x is ', num2str(phi_x1 * phi_x), ...
				'; phi_x1 is ', num2str(phi_x1), '; phi_x is ', num2str(phi_x)]);
        end
        K_approx(i, j) = k_approx;
        % True SE kernel
        dist_sq = sum((x1 - x).^2);
        K_true(i, j) = sigma2 * exp(-dist_sq / (2 * lengthScale^2));
    end
end

% Plot the approximate SE kernel
figure;
surf(X, Y, K_approx);
title('Approximated SE Kernel');
xlabel('x1');
ylabel('x2');
zlabel('Kernel Value');
shading interp;
colorbar;

% Plot the true SE kernel
figure;
surf(X, Y, K_true);
title('True SE Kernel');
xlabel('x1');
ylabel('x2');
zlabel('Kernel Value');
shading interp;
colorbar;

% Plot the difference between the approximate and true kernels
figure;
surf(X, Y, abs(K_approx - K_true));
title('Difference between Approximate and True SE Kernels');
xlabel('x1');
ylabel('x2');
zlabel('Difference');
shading interp;
colorbar;

function [eigenval,eigenfun,eigenfun_dx,NN] = domain_cartesian_dx(m,d,L)
%% domain_cartesian_dx - Laplace operator eigendecomposition in a hypercube
% 
% Syntax:
%  [eigenfun,eigenval,NN] = domain_cartesian(m,d,L)
%
% In:
%   m  - Number of eigenfunctions
%   d  - Dimensionality
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   eigenval    - Function handle: eigenval(n)
%   eigenfun    - Function handle: eigenfun(n,x)
%   eigenfun_dx - Function handle: eigenfun_dx(n,x)
%   NN          - Indices to evaluate the handles at
% 
% Description:
%   This code returns the eigendecomposition of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
% Copyright (c) 2014-2023 Arno Solin

  % If domain is not centered, center it
  if size(L,1)>1
    L = (max(L,[],1)-min(L,[],1))/2;
  end

  % This is stupid, but at least we should get enough 
  % of basis function to choose from
  N = ceil(m^(1/d)*L/min(L)); % BB: N is number of basis function, seems to be 512
  
  % Combined eigenfunction indices (checked the numbers)
  NN = ndgridm(N);

  % Define eigenvalues of the negative Laplacian in ND 
  % s.t. Dirichlet boundary. This forms an orthonormal basis.
  eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
    
  % Sort and take only m most important eigenfunctions
  [~,ind] = sort(eigenval(NN)); 
  disp(['ind length is ',num2str(length(ind))]); 
  disp(['m is ', num2str(m)]);
  NN = NN(ind(1:m),:);  

  % Define eigenfunction of the negative Laplacian in ND 
  % s.t. Dirichlet boundary. This forms an orthonormal basis.
  eigenfun = @(n,x) laplace_eig_cart_dirichlet(n,x,L); % n  - Eigenfunction indices;  x  - Spatial locations [x1 x2]

  % Define derivative of the eigenfunction of the negative 
  % Laplacian in ND s.t. Dirichlet boundary 
  eigenfun_dx = @(n,x,di) laplace_eig_cart_dirichlet_dx(n,x,di,L);
  
  
end


function [v]=laplace_eig_cart_dirichlet(n,x,L)
%% laplace_eig_cart_dirichlet - Laplace operator eigenfunctions in a hypercube
% 
% Syntax:
%  [v] = laplace_eig_cart_dirichlet(n,x,L)
%
% In:
%   n  - Eigenfunction indices
%   x  - Spatial locations [x1 x2]
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   v - The evaluated value
% 
% Description:
%   This code calculates the eigenvectors of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
%   The corresponding eigenvalues can be calculated by
% 
%     eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
%
% Copyright (C) 2012 Arno Solin
%

  % Allocate space
  v = ones(size(x,1),size(n,1));

  % Evaluate eigenfunctions
  for j=1:size(n,2)
    for i=1:size(n,1)
      v(:,i) = v(:,i) .* 1./sqrt(L(j)) .* ...
        sin(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
    end
  end
  
  
%   % Evaluate eigenfunctions
%   if size(x,2)==1
%       for i=1:numel(n)
%           v(:,i) = sqrt(1./L)*sin(pi*n(i)*(x(:)+L)/2/L);
%       end
%   else
%       for i=1:size(n,1)
%           % Eigenfunctions for x in Omega and n = (n1,n2,...nn), ni = 1,2,...,Ni
%           v(:,i) = prod(bsxfun(@times,sqrt(1./L), ...
%               sin(pi*bsxfun(@times,n(i,:)./L,bsxfun(@plus,x,L))/2)),2);
%           if all(n(i,:)==0)
%               v(:,i) = ones(size(x,1),1);
%           end
%       end
%   end

end

function [v]=laplace_eig_cart_dirichlet_dx(n,x,di,L)
%% laplace_eig_cart_dirichlet_dx - Derivative of Laplace operator eigenfunctions in a hypercube
% 
% Syntax:
%  [v] = laplace_eig_cart_dirichlet_dx(n,x,di,L)
%
% In:
%   n  - Eigenfunction indices
%   x  - Spatial locations [x1 x2]
%   di - Differentiate w.r.t. this dimension d/dx_i
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   v - The evaluated value
% 
% Description:
%   This code calculates the eigenvectors of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
%   The corresponding eigenvalues can be calculated by
% 
%     eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
%
% Copyright (C) 2015 Arno Solin
%

  % Allocate space
  v = ones(size(x,1),size(n,1));

  % Evaluate eigenfunctions
  for j=1:size(n,2)
    if j==di % Differentiate w.r.t this dimension
      for i=1:size(n,1)
        %v(:,i) = v(:,i) .* 1./sqrt(L(j)) .* ...
        %  cos(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
        
        % This is the actual derivative  
        v(:,i) = v(:,i) .* pi*n(i,j)/(2*L(j)*sqrt(L(j))) .* ...
          cos(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
      
        %v(:,i) = v(:,i).*prod(bsxfun(@times,sqrt(1./L(j)), ...
        %   cos(pi*bsxfun(@times,n(i,j)./L(j),bsxfun(@plus,x(:,j),L(j)))/2)),2);
      end
    else % Do not differentiate w.r.t this dimension
      for i=1:size(n,1)
        v(:,i) = v(:,i) .* 1./sqrt(L(j)) .* ...
          sin(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
        %v(:,i) = v(:,i).*prod(bsxfun(@times,sqrt(1./L(j)), ...
        %   sin(pi*bsxfun(@times,n(i,j)./L(j),bsxfun(@plus,x(:,j),L(j)))/2)),2);       
        %if all(n(i,j)==0)
        %  %v(:,i) = ones(size(x,1),1);
        %end
      end
    end
  end
  
end

function NN = ndgridm(N)
%% ndgridm - Expand index hypercude
%
% Syntax:
%  [NN] = ndgridm(N)
%
% In:
%   N  - Vector of max indices
%      
% Out:
%   NN - Matrix of index combinations
%
% Description:
%   A more felxible variant of 'ndgrid'. This functions gives combinations
%   of the indices in N such that, for example, for a 2D case we get
%   (1,1),(1,2),...,(1,N2),(2,1),...,(2,N2),(...),(N1,1),...,(N1,N2).
%   This function works for any dimension.
%
% Copyright (C) 2014 Arno Solin
%

  % Allocate space for indices
  NN = zeros(prod(N),numel(N));

  % For each level/diemsion
  if numel(N)==1

     % The lowest level
     NN(:,1) = (1:N)';
     
  else

    % This level
    n = 1:N(1);

    % Recursive call
    nn = ndgridm(N(2:end));

    % Assign values
    NN(:,1)     = kron(n,ones(1,prod(N(2:end))))';
    NN(:,2:end) = repmat(nn,[N(1) 1]);
   
  end

end
