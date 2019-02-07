%% function that returns the negative log likelihood of a logistic regression model given parameter values and data

function [out] = negloglikGLM (paras,response,temperature,dummy1,dummy2)

%% 'paras' is the vector of model parameters.
%% we want to estimate values for 'paras' that MINIMISE the value of the function 'negloglik'

% first to fourth entries of 'paras' vector are beta_0 to beta_3
% 'dummy1' and 'dummy2' are dummy variables for polymer grade

% work out the linear predictor:
linear_predictor = paras(1) + paras(2).*temperature + paras(3).*dummy1 + paras(4).*dummy2; % here we work out the predicted mean of the response Y based on our model

% send the linear predictor through the logistic function to get
% probabilities:
% TIP: P(failure) = 1-P(success)
probs1 = ones(1,length(response))./(1+exp(-linear_predictor));
probs = response.*probs1 + (1-response).*(1-probs1);
out = sum(-log(probs)); % we return the negative log-likelihood 
%% This is convenient, as products turn into sums when we take logs. This means we do not have to worry about 'underflow errors'.
%% these occur when a number becomes smaller than the computer can handle, which happens frequently when multiplying many small probability densities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%