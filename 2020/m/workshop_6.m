%% Workshop 6
%% 1 (DRILLROCK)
delimiterIn = ',';
A = importdata('DRILLROCK.csv',delimiterIn);
% how many data points are there?
A.data
size(A.data)
% plot scatterplots
clf
x = A.data(:,1); % depth data
y = A.data(:,2); % time data
plot(x,y,'.')
title('Scatterplot') % title for plot
xlabel('depth') % x-axis label
ylabel('time') % y-axis label

% 1.2
Y = y;
% define the matrix X. Don't forget about the intercept.
X = [ones(17,1) x];
betahat = inv(X'*X)*(X'*Y);
sigmahat = sqrt((Y-X*betahat)'*(Y-X*betahat)/(17-2));

% define values for parameters that we use at the start of the optimisation.
% We have 3 parameters: sigma, beta_0, beta_1
% x0 is the same length as the vector 'paras' in the function 'negloglik'.
x0 = [sqrt(var(y)) 1 0];
% check that computing the negative log-likelihood using 'negloglik.m' works:
negloglik([1 2 0],1,2) % uses sigma=1, beta_0=2 and beta_1=0 for data: Y=1 and x=2.
negloglik(x0,y,x) % uses starting values of parameters and all data.
% maximum likelihood fit:
fun = @(s)negloglik(s,y,x); % this sets up the function we want to minimise.
s = fminsearch(fun,x0); % here we performe the optimisation
% (i.e. find parameters for which negloglik is minimal)
% The fitted parameters are stored in 's':

% create a table for the data (x,y, defined above):
data = table(x,y,'VariableNames',{'depth','time'});
% show first few rows of table
data(1:5,:)
% fit simple linear model
m1 = fitlm(data, 'time~depth');

% 1.3
% use parameter estimates and data to find residuals.
resids = Y-X*betahat;
% plot a histogram of the residuals
clf
subplot(1,2,1)
hist(resids)
xlabel('Residuals')
ylabel('Frequency')
% plot the residuals versus the fitted values
fitted = X*betahat;
subplot(1,2,2)
plot(fitted,resids,'.')
xlabel('Fitted values') % x-axis label
ylabel('Residuals') % y-axis label

% plot distribution of residuals (for outliers etc.)
clf
subplot(2,2,1)
plotResiduals(m1)
% plot to check normality
subplot(2,2,2)
plotResiduals(m1,'probability')
% residuals versus fitted values (check for homoscedasticity)
subplot(2,2,3)
plotResiduals(m1,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m1,'lagged') % want no trend in this!

% 1.4. Check the model fit visually
% check model fit:
clf
% plot the raw data:
plot(x,y,'.')
hold on
% plot model fit line using the estimates parameters:
xx = [0 400];
yy = [m1.Coefficients{1,1}+0*m1.Coefficients{2,1} m1.Coefficients{1,1}+400*m1.Coefficients{2,1}];
plot(xx,yy,'k','LineWidth',2)
xlabel('depth') % x-axis label
ylabel('time') % y-axis label
title('Model fit line')

% 1.5 
% re-calculate the error estimate
sigmahat = sqrt((Y-X*betahat)'*(Y-X*betahat)/(17-2));
% calculate the quantity SS_xx
ssxx = sum((x-mean(x)).*(x-mean(x)));
% calculate the T-value
T = (betahat(2)-0)/(sigmahat/sqrt(ssxx));
% notice how you can also obtain this from the first
% and second columns of the Matlab output:
m1.Coefficients{2,1}/m1.Coefficients{2,2};
% obtain the p-value (two-tailed):
P = 2*(1-tcdf(T,17-2));