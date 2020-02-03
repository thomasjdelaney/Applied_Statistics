%% lab02_lab.m
%% Parameter Estimation
%% 1. Maximum Likelihood
%% Bernoulli distribution
p_true = 0.75;
pd = makedist('Binomial', 'p', p_true);
x = random(pd, 1, 10);

% define the likelihood for the binomial distribution
likelihood = @(p) prod(p.^x .* (1-p).^(1-x));
% define the log-likelihood for the binomial distribution
loglikelihood = @(p) sum(x .* log(p) + (1-x) .* log(1-p));

% compute likelihood
p = 0: 0.001: 1;
likelihoods = zeros(1, length(p));
loglikelihoods = zeros(1, length(p));
for i = 1:length(p)
    likelihoods(i) = likelihood(p(i));
    loglikelihoods(i) = loglikelihood(p(i));
end

% plot
clf;
yyaxis left
plot(p, likelihoods, '-', 'Linewidth', 2)
ylabel('L') % y-axis label
yyaxis right
plot(p, loglikelihoods, '-', 'Linewidth', 2)
hold on

% Find p that maximises L (i.e., minimises -L)
negL = @(x) likelihood(x) * -1;
p_L = fminsearch(negL, 0); % initial value of p is 0

% Find p that maximises logL (i.e., minimises -logL)
neglogL = @(x) loglikelihood(x) * -1;
% fminsearch(neglogL, 0)
p_logL = fminsearch(neglogL, 0);

p_hat = mean(x);
stem(p_hat, 1, '-y', 'Linewidth', 5);
hold on
stem(p_L, 1, '-.g', 'Linewidth', 3);
hold on
stem(p_logL, 1, '--r', 'Linewidth', 3);
hold on
stem(p_true, 1, '-.', 'Linewidth', 2);
legend('Likelihood', 'logLikelihood', '$p_{MLE}$', '$p_{L}$', '$p_{logL}$', '$p_{true}$', 'Location','northwest', 'Interpreter', 'latex')
title('MLE for Binomial') % title for plot
xlabel('p') % x-axis label
ylabel('logL') % y-axis label

% Matlab's MLE function
mle(x,'distribution','Bernoulli')

%% Normal distribution
mu_true = 0;
sigma_true = 1;
pd = makedist('Normal', 'mu', mu_true, 'sigma', sigma_true);
x = random(pd, 1, 10);

% define the likelihood for the binomial distribution
likelihood = @(mu, sigma) prod(normpdf(x, mu, sigma));
% define the log-likelihood for the binomial distribution
loglikelihood = @(mu, sigma) sum(log(normpdf(x, mu, sigma)));
% plot
clf;
fcontour(likelihood, [-1 1 0 2],'Fill','on','MeshDensity',200);

hold on
plot(mu_true, sigma_true, 'dr', 'MarkerSize', 10, 'MarkerFaceColor',"r")
hold on
% Find p that maximises L (i.e., minimises -L)
negL = @(ms) likelihood(ms(1), ms(2)) * -1;
msL = fminsearch(negL, [0.5, 0.5]);

plot(msL(1), msL(2), 'sg', 'MarkerSize', 15, 'MarkerFaceColor',"g")
hold on
% Find p that maximises logL (i.e., minimises -logL)
neglogL = @(ms) loglikelihood(ms(1), ms(2)) * -1;
% fminsearch(neglogL, 0)
mslL = fminsearch(neglogL, [0.5, 0.5]);

plot(mslL(1), mslL(2), 'om', 'MarkerSize', 8, 'MarkerFaceColor',"m")
hold on
% MLE estimate
plot(mean(x), std(x, 1), '*k', 'MarkerSize', 10, 'MarkerFaceColor',"k")
title('Likelihood of Normal distribution') % title for plot
xlabel('$\mu$', 'Interpreter', 'latex') % x-axis label
ylabel('$\sigma$', 'Interpreter', 'latex') % y-axis label

% Matlab's MLE function
mle(x,'distribution','Normal')

%% Compute MLE by hand
theta = 3.5;

n = 10000;
x = rand(1, n) .^ (-1/theta);
n / sum(log(x));

%% 2. Methods of Moments
%% MoM estimate for parameters of Poisson distribution
% generate Poisson random sample
lambda = 10;
n = 10000;
x = poissrnd(lambda, 1, n);
% MoM estimate of parameter lambda
mean(x);

% increasing sample size
ns = round(logspace(0, 4, 20));
s = 1000;
% compute the standard deviation of the estimates for different sample sizes
stdevs = zeros(1, length(ns));
for i = 1:length(ns)
    n = ns(i);
    estimates = zeros(1, s);
    for j = 1:s
        x = poissrnd(lambda, 1, n);
        estimates(j) = mean(x);
    end
    stdevs(i) = std(estimates);
end
% plot standard deviation (spread) of estimates against the sample size
clf;
plot(ns, stdevs, 'o');
hold on
plot(ns, ns.^(-0.5), '--');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('$n$', 'Interpreter', 'latex') % x-axis label
ylabel('$\sigma_{\hat{\lambda}}(n)$', 'Interpreter', 'latex') % y-axis label

%% 3. Unbiased Estimators
% true value
sigma = 5;

% sample
n = 10000;
x = exprnd(sigma, 1, n) .* (2 * binornd(1, 0.5, 1, n) - 1);
% pdf
clf;
histogram(x, 'Normalization',"pdf");
hold on
xx = min(x):0.1:max(x);
plot(xx, exp(- abs(xx) / sigma) / 2 / sigma, '-r', 'LineWidth', 2)

close all
% estimators
sigma_mle = @(x) sum(abs(x)) / length(x);
sigma_mom = @(x) (sum(x.^2) / 2 / length(x))^0.5;
% estimates distributions
S = 20000;
ns = [5 10, 100];

estimates_mle = zeros(length(ns), S);
estimates_mom = zeros(length(ns), S);
for i = 1:length(ns)
    n = ns(i);
    for j = 1:S
        x = exprnd(sigma, 1, n) .* (2 * binornd(1, 0.5, 1, n) - 1);
        estimates_mle(i, j) = sigma_mle(x);
        estimates_mom(i, j) = sigma_mom(x);
    end
end

% set(gcf, 'Position', [0 0 1200 300])
clf;
for i = 1:length(ns)
    figure(i)
    histogram(estimates_mle(i, :), 'Normalization',"pdf");
    hold on
    histogram(estimates_mom(i, :), 'Normalization',"pdf");
    hold on;
    stem(sigma, 1, '--r', 'Linewidth', 4);
    hold on
    stem(mean(estimates_mle(i, :)), 1, '--', 'Linewidth', 2);
    hold on
    stem(mean(estimates_mom(i, :)), 1, '-.y', 'Linewidth', 2);
    legend('mle', 'mom', 'true', 'mean mle', 'mean mom');
    title(sprintf('sample size n = %d', ns(i)))%, 'Interpreter', 'latex')
    xlabel('Estimate');
    ylabel('PDF');
    % xlim([0 4]); % uncomment this if you want all plots on the same scale
end

%% 5. Estimate the heights of the tallest students
true_mean = 167;
true_std = 10;
% True mean of the top 10% of the distribution of heights
m = true_mean;
v = true_std^2;
th_q09 = norminv(0.9, m, v^0.5);

syms y(h)
y(h) = 10 * h * exp(-(h - m)^2 / 2 / v) / (2 * pi * v)^0.5;

% sample the heights of n students
n = 10;
x = normrnd(true_mean, true_std, 1, n);
% h_s first estimator: the mean height of the top 10% highest student in the sample
% compute the empirical quantile 0.9
q09 = quantile(x, 0.9);
h_true = vpaintegral(y, th_q09, th_q09 * 10);

tallest = x(x > q09);
hs = mean(tallest);

% h_N second estimator: the mean of the top 10% of the distribution of
% heights, whose parameters are estimated from the sample
m = mean(x)
v = std(x)^2
th_q09 = norminv(0.9, m, v^0.5)
syms y(h)
y(h) = 10 * h * exp(-(h - m)^2 / 2 / v) / (2 * pi * v)^0.5;
hN = vpaintegral(y, th_q09, th_q09 * 10)

% estimates distributions
S = 100; % fast, but not enough samples
% S = 1000; % ...this is quite slow!
ns = [10 20, 100];
estimates_hs = zeros(length(ns), S);
estimates_hN = zeros(length(ns), S);
for i = 1:length(ns)
    n = ns(i);
    for j = 1:S
        x = normrnd(true_mean, true_std, 1, n);
        
        % hs
        q09 = quantile(x, 0.9);
        tallest = x(x > q09);
        hs = mean(tallest);
        % hN
        m = mean(x);
        v = std(x)^2;
        th_q09 = norminv(0.9, m, v^0.5);
        syms y(h)
        y(h) = 10 * h * exp(-(h - m)^2 / 2 / v) / (2 * pi * v)^0.5;
        hN = vpaintegral(y, th_q09, th_q09 * 10);
        estimates_hs(i, j) = hs;
        estimates_hN(i, j) = hN;
    end
end

close all
% set(gcf, 'Position', [0 0 1200 300])
clf;
for i = 1:length(ns)
    figure(i)
    histogram(estimates_hs(i, :), 'Normalization',"pdf");
    hold on
    histogram(estimates_hN(i, :), 'Normalization',"pdf");
    hold on
    stem(184.55, 0.3, '--r', 'Linewidth', 4); % h_true doesn't work..?
    hold on
    stem(mean(estimates_hs(i, :)), 0.3, '--', 'Linewidth', 2);
    hold on
    stem(mean(estimates_hN(i, :)), 0.3, '-.y', 'Linewidth', 2);
    legend('PDF of $\hat{h}_s$', 'PDF of $\hat{h}_N$', 'true $h$', 'mean $\hat{h}_s$', 'mean $\hat{h}_N$', 'Interpreter', 'latex');
    title(sprintf('sample size n = %d', ns(i)))%, 'Interpreter', 'latex')
    xlabel('Estimate');
    ylabel('PDF');
end