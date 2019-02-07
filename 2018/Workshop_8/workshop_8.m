%% Workshop 8
%% Question 1.0 Different types of preditctors
% One solution:
x = random('Uniform',0,10,[1000,1]); % predictor
beta0 = 2.3; % intercept
beta1 = 1.5; % effect of predictor
epsilon = random('Normal',0,2.1,[1000,1]); % error vector N(0,2.1)
% define response 
y = beta0 + beta1.*x + epsilon;

% check the scatterplot of these data
clf
plot(x,y,'.'); xlabel('x'); ylabel('y');

% fit a linear model
data = table(x, y, 'VariableNames', {'predictor', 'response'});
m1 = fitlm(data, 'response~predictor');
% plot distribution of residuals
clf
subplot(2,2,1); plotResiduals(m1); % histogram
subplot(2,2,2); plotResiduals(m1, 'probability'); % Q-Q to check normality
subplot(2,2,3); plotResiduals(m1, 'fitted'); % residuals vs fitted values
subplot(2,2,4); plotResiduals(m1, 'lagged'); % auto correlation via lagged residuals

%% Question 1.2
% qualitative predictor with 3 levels
x1 = random('Uniform', 0, 10, [1000, 1]); % quantitative predictor
x2 = zeros(1000, 1); % qualitative predictor
x2(1:400) = 0; % level 0
x2(401:700) = 1; % level 1
x2(710:1000) = 2; % level 2
beta0 = 2.3; % intercept
beta1 = 1.5; % effect of predictor x1
beta2 = 6; % effect of level 1 of x2
beta3 = -5.5; % effect level 2 of x2
epsilon = random('Normal', 0, 3.5, [1000, 1]); % error vector N(0, 3.5)
% create some dummy variables for simulating the data
dum1 = zeros(1000, 1);
dum2 = zeros(1000, 1);
dum1(x2==1) = 1;
dum2(x2==2) = 1;
y = beta0 + beta1.*x1 + beta2.*dum1 + beta3.*dum2 + epsilon;

% scatterplots
clf
subplot(2,1,1)
plot(x1,y,'.'); xlabel('x1'); ylabel('y');
subplot(2,1,2)
plot(x2,y,'.'); xlabel('x2'); ylabel('y');

% fitting a linear model
data = table(x1, x2, y, 'VariableNames', {'predictor1', 'predictor2', 'response'});
m1 = fitlm(data, 'response~predictor1+predictor2');

% plot residuals
clf
subplot(2,1,1)
plot(x1, m1.Residuals{:,1}, '.'); xlabel('x1'); ylabel('raw residuals');
subplot(2,1,2)
plot(x2, m1.Residuals{:,1}, '.'); xlabel('x2'); ylabel('raw residuals');

% tell matlab that x2 is a qualitative predictor
data.predictor2 = nominal(data.predictor2);

% fit the model again
m2 = fitlm(data, 'response~predictor1+predictor2');

% plot residuals again
clf
subplot(2,1,1)
plot(x1, m2.Residuals{:,1}, '.'); xlabel('x1'); ylabel('raw residuals');
subplot(2,1,2)
plot(x2, m2.Residuals{:,1}, '.'); xlabel('x2'); ylabel('raw residuals');

%% Question 1.3 Interaction Terms
% One solution:
 
x1 = random('Uniform',0,10,[1000,1]); % quantitative predictor
x2 = zeros(1000,1); % this will be our qualitative predictor
x2(1:500) = 0; % level 0
x2(501:1000) = 1; % level 1
beta0 = 2.3; % intercept
beta1 = 1.5; % baseline effect of predictor x1
beta2 = 6; % baseline effect of level 1 of x2
beta3 = -3.5; % interaction between x1 and level 1 of x2
epsilon = random('Normal',0,5,[1000,1]); % error vector N(0,5)
% construct response:
y = beta0 + beta1.*x1 + beta2.*x2 + beta3.*x1.*x2 + epsilon;

% scatter plot
clf
plot(x1,y,'.'); xlabel('x1'); ylabel('y')

% fit a linear model
data = table(x1,x2,y,'VariableNames',{'predictor1','predictor2','response'});
m1 = fitlm(data, 'response~predictor1+predictor2');

% plot residuals
clf
plotResiduals(m1,'fitted')

% fit the model again with a qualitative predictor and an interaction term
data.predictor2 = nominal(data.predictor2);
m2 = fitlm(data, 'response~predictor1*predictor2');

% residual plots again
clf
subplot(2,2,1); plotResiduals(m2)
subplot(2,2,2); plotResiduals(m2,'probability') % Q-Q plot to check normality
subplot(2,2,3); plotResiduals(m2,'fitted') % residuals versus fitted values
subplot(2,2,4); plotResiduals(m2,'lagged') % auto-correlation (via lagged residuals)


%% Question 1.4 Polynomials of Predictors
% One solution:
 
x = random('Uniform',-10,10,[1000,1]); % predictor
beta0 = 2.3; % intercept
beta1 = 1.5; % effect of predictor
beta2 = 0.3; % effect of square of predictor
epsilon = random('Normal',0,4.3,[1000,1]); % error vector N(0,4.3)
y = beta0 + beta1.*x + beta2.*x.*x + epsilon;

% scatterplot
clf
plot(x, y, '.'); xlabel('x'); ylabel('y');

% fit a linear model
data = table(x,y,'VariableNames',{'predictor','response'});
m1 = fitlm(data, 'response~predictor');

% plot the residuals
clf
subplot(2,2,1)
plotResiduals(m1)
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m1,'probability')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m1,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m1,'lagged')

% fit a model with a quadratic term
data2 = table(x,x.^2,y,'VariableNames',{'x','xsquared','response'});
m2 = fitlm(data2, 'response~x+xsquared');

% plot the residuals again
clf
subplot(2,2,1)
plotResiduals(m2)
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m2,'probability')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m2,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m2,'lagged')

%% skipping to Question 2: Common Pitfalls
%% Question 2.2 Violated model assumption
% 2.2.1 Non-normal error/noise distribution
x = random('Uniform',0,10,[1000,1]); % predictor
beta0 = 2.3; % intercept
beta1 = 1.5; % effect of predictor
epsilon = random('Exponential',5,[1000,1]); % error vector
y = beta0 + beta1.*x + epsilon;

% scatterplot
clf
plot(x, y, '.'); xlabel('x'); ylabel('y');

% fit a linear model
data = table(x,y,'VariableNames',{'predictor','response'});
m1 = fitlm(data, 'response~predictor');

% residual plots
clf
subplot(2,2,1)
plotResiduals(m1)
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m1,'probability')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m1,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m1,'lagged')

% 2.2.2 Non-constant variance of predictors
x = random('Uniform',0,10,[1000,1]); % predictor
beta0 = 2.3; % intercept
beta1 = 1.5; % effect of predictor
epsilon = random('Normal',0,0.5*x,[1000,1]); % variance is proportional to predictor
y = beta0 + beta1.*x + epsilon;

% scatterplot
clf
plot(x, y, '.'); xlabel('x'); ylabel('y');

% fit a linear model
data = table(x,y,'VariableNames',{'predictor','response'});
m1 = fitlm(data, 'response~predictor');

% plot distribution of residuals (for outliers etc.)
clf
subplot(2,2,1)
plotResiduals(m1)
% Q-Q plot to check normality
subplot(2,2,2)
plotResiduals(m1,'probability')
% residuals versus fitted values
subplot(2,2,3)
plotResiduals(m1,'fitted')
% auto-correlation (via lagged residuals)
subplot(2,2,4)
plotResiduals(m1,'lagged')