%% Revision Workshop
%% Question 1: Outliers
x = load('outliers.txt');
x = sort(x(:));
[y, idx, outl] = deleteoutliers(x, 0.2); % delete outliers
if isempty(idx)
    disp('No outliers.')
else
    fprintf('Num outliers = %d\n', length(idx))
end

%% Question 2: Polynomial regression
xy = load('regr.txt');
x = xy(:,1); y = xy(:,2);
n = size(xy, 1);
fprintf('Number of paired observations n = %d. \n', n);

sg = 0.01; % significance level
kmax = 20; % max degree for the polynomial is set to 20
k = [0:kmax];
X = repmat(x,1,kmax+1).^repmat(k,n,1); % form matrix of functions 1, x^2
X(:,1) = X(:,1)/norm(X(:,1)); % normalise first column
for ki = 2:kmax+1 % orthonormalise columns
    Sum = 0;
    for kk = 1:ki-1
        Sum = Sum + (X(:,ki)'*X(:,kk))*X(:,kk);
    end
    X(:,ki) = X(:,ki) - Sum;
    X(:,ki) = X(:,ki)/norm(X(:,ki));
end
[b,bint]=regress(y,X,sg); % apply regression model

        