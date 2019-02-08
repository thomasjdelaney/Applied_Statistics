%% 1
% a
est=xlsread('aircraft_cond.xlsx');% read xlsx file
est=est(:); % convert matrix to vector
id1 = find(isnan(est)); % find empty elements in the vector
est(id1) = []; % delete empty elements - now vector contains 199 non-empty values
[muhat,muci] = expfit(est); % point estimate for mean (muhat) and confidence
[muhat,muci] = expfit(est, 0.01); % estimates with 99% confidence interval

% b
[phat,pci] = mle('exp',est);

% c
% In lecture, lambda (rate) was estimated

% d 
lambda_bar = 1/muhat;

% e and f
% formlae for standard error, and z-score coefficient available in lecture

%% 2
% a
t = 0:0.1:10;
h_e = exppdf(t,0.7)./(1-expcdf(t,0.7));
h_w = wblpdf(t,2,0.5)./(1-wblcdf(t,2,0.5));
plot(t,h_e,'*');
hold on
plot(t,h_w,'-');

% c
failure_time = [1,2,3,4,5,6];
freq = [6, 3, 2, 2, 2, 1];
[params_est, params_ci] = wblfit(freq);