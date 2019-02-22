%% For workshop number 4 of the applied statistics course
%% Question 1
%%
drugtest = xlsread('DRUGTEST.XLS');
% H0 is that the death rates are the same, HA is that they are different
suff = drugtest(:,1);
clin = drugtest(:,2);
clin(isnan(clin)) = []; % remove nans
sum_suff = sum(suff); sum_clin = sum(clin);
n_suff = length(suff); n_clin = length(clin);

% null hypothesis
all_data = [suff; clin];
n = length(all_data);
sum_all_data = sum(all_data);
lambda_mle_0 = n/sum_all_data;
log_like_lambda_mle_0 = n*log(lambda_mle_0) - sum_all_data*lambda_mle_0;

% alternative hypothesis
lambda_mle_suff = n_suff/sum_suff;
lambda_mle_clin = n_clin/sum_clin;
log_like_lambda_mle_suff = n_suff*log(lambda_mle_suff) - sum_suff*lambda_mle_suff;
log_like_lambda_mle_clin = n_clin*log(lambda_mle_clin) - sum_clin*lambda_mle_clin;
log_like_lambda_mle_A = log_like_lambda_mle_suff + log_like_lambda_mle_clin;

% ratio test
[h,pValue,stat] = lratiotest(log_like_lambda_mle_A,log_like_lambda_mle_0,1);
% ratio test by hand
chi_squared = 2*(log_like_lambda_mle_A - log_like_lambda_mle_0);
p_value = 1 - chi2cdf(chi_squared, 1);
% There's one degree of freedom because the alternative hypothesis has
% two parameters, whereas the null hypothesis only has one
% The low p-value implies that the null hypothesis does not hold.

%% Question 2 (a)
%%
% Need to convert SILICA.XLS to a csv because SILICA.XLS contains unicode
% characters
silica = csvread('SILICA.csv', 1);
% It's a two tailed test because we are looking for any effect.
% H0 both sample sets have the same mean
% HA both sample sets have different means
without = silica(:,1); with = silica(:,2);
[h, p, ci, stats] = ttest2(without, with, 'Alpha', 0.1); % rejected!

% z-stat
n = length(without);
mean_without = mean(without); mean_with = mean(with);
var_without = var(without); var_with = var(with);
z_stat = (mean_without - mean_with)/sqrt((var_without + var_with)/n);
rej_region = norminv(0.95, 0, std([without; with]));

% same result, different p value, a lot of samples makes the z-test
% reasonable

%% Question 2 (b)
%%
% H0: Flexed has mean = mean of unflexed
% HA: Flexed has mean <= mean of unflexed
% I suppose it's a one tailed test
n=11;
mean_flex = 59; mean_unflex = 43;
std_flex = 4; std_unflex = 2;

% t(n1 + n2 - 2)
rejection_region = tinv(0.95, n+n-2);

mean_diff = mean_flex - mean_unflex;
full_var = (10*(std_flex*std_flex) + 10*(std_unflex*std_unflex))/20;
t_val = mean_diff/sqrt(full_var*(2/11));
% t_val > rejection_region, therefore the null hypothesis is rejected
% i.e. mean flexed > mean unflexed

std_flex = 10; std_unflex = 15;
mean_diff = mean_flex - mean_unflex;
full_var = (10*(std_flex*std_flex) + 10*(std_unflex*std_unflex))/20;
t_val = mean_diff/sqrt(full_var*(2/11));

%% Question 3 (a)
%%
shallow = csvread('SHALLOW.csv', 1);
actual = shallow(:,1);
predicted = shallow(:,2);
[h, p, ci, stats] = ttest(actual, predicted, 'alpha', 0.05);

%% Qestion 4 (a)
%%
% H0: mean client-flow = mean server-flow
% HA: mean client-flow > Mean server-flow
n_server = 40; vul_server = 20;
n_client = 54; vul_client = 41;
n_total = n_server + n_client;
vul_total = vul_server + vul_client;

phat_server = vul_server/n_server;
phat_client = vul_client/n_client;
phat_total = vul_total/n_total;
var_total = phat_total*(1 - phat_total);

z_stat= (phat_server - phat_client)/sqrt(var_total*(1/n_server + 1/n_client));
rejection_region = norminv(0.01, 0, sqrt(var_total));