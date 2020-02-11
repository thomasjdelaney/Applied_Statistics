%% lab03_lab.m
%% Confidence Intervals
%% 1. z-CI and t-CI for the mean of a Normal RV

% Normal RV
true_mu = 0;
true_sigma = 1;
pd = makedist('Normal', 'mu', true_mu, 'sigma', true_sigma);

% sampling distribution
% generate samples and estimates
ns = [5 20 100]; % size of samples
S = 10000; % number of samples
estimates = zeros(length(ns), S);
sample_stdevs = zeros(length(ns), S);
for i = 1:length(ns)
    n = ns(i);
    for j = 1:S
        % sample from the normal RV
        x = random(pd, 1, n);
        % compute the estimate for the sample
        estimates(i, j) = mean(x);
        sample_stdevs(i, j) = std(x);
    end
end

% plot the sampling distributions
clf;
for i = 1:length(ns)
    x = estimates(i, :);
    histogram(x, 'Normalization', "pdf", 'EdgeColor','none');
    hold on
end
stem(true_mu, 5, '--k', 'Linewidth', 2);
hold on
xlim([-1.5 1.5]);
title('Sampling distribution for the mean of a Normal RV') %, 'Interpreter', 'latex')
legend('sample size n = 5', 'sample size n = 20', 'sample size n = 100', 'true value', 'FontSize', 10.0);
xlabel('Estimate $\bar{X}$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);

% compute empirical Confidence Intervals
clf
alpha = 0.05;
eCIs = zeros(length(ns), 2);
for i = 1:length(ns)
    eCIs(i, 1) = quantile(estimates(i, :), alpha / 2);
    eCIs(i, 2) = quantile(estimates(i, :), 1 - alpha / 2);
    % lengths
    eCIs(i, 2) - eCIs(i, 1)
    % plot pdf
    x = estimates(i, :);
    histogram(x, 'Normalization', "pdf", 'EdgeColor','none');
    hold on
end

stem(true_mu, 5, '--k', 'Linewidth', 2);
hold on
stem(eCIs(1, 1), 5, '--b', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(1, 2), 5, '--b', 'Linewidth', 3,'HandleVisibility','off');
hold on
stem(eCIs(2, 1), 5, '--r', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(2, 2), 5, '--r', 'Linewidth', 3,'HandleVisibility','off');
hold on
stem(eCIs(3, 1), 5, '--y', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(3, 2), 5, '--y', 'Linewidth', 3,'HandleVisibility','off');
xlim([-1.5 1.5]);

% compute the z-CI (Lecture 3, Eq. 12) for each sample
zCIs1 = zeros(length(ns), S); % lower edges
zCIs2 = zeros(length(ns), S); % upper edges
for i = 1:length(ns)
    zCIs1(i, :) = estimates(i, :) - norminv(1 - alpha / 2) * 1 / (ns(i))^0.5;
    zCIs2(i, :) = estimates(i, :) - norminv(alpha / 2) * 1 / (ns(i))^0.5;
    % average lengths
    mean(zCIs2(i, :) - zCIs1(i, :), 2)
end

% fraction of intervals that contain the true mean, for the various sizes n
mean(zCIs1 < true_mu & true_mu < zCIs2, 2) % mean over the second dimension (axis)

% compute the t-CI (Lecture 3, Eq. 13) for each sample
tCIs1 = zeros(length(ns), S); % lower edges
tCIs2 = zeros(length(ns), S); % upper edges
for i = 1:length(ns)
    tCIs1(i, :) = estimates(i, :) - tinv(1 - alpha / 2, ns(i) - 1) * sample_stdevs(i, :) / (ns(i))^0.5;
    tCIs2(i, :) = estimates(i, :) - tinv(alpha / 2, ns(i) - 1) * sample_stdevs(i, :) / (ns(i))^0.5;
    % average sizes
    mean(tCIs2(i, :) - tCIs1(i, :), 2)
end

% fraction of intervals that contain the true mean, for the various sizes n
mean(tCIs1 < true_mu & true_mu < tCIs2, 2) % mean over the second dimension (axis)

%% 2. CI for the mean of a Bernoulli RV & the Central Limit Theorem

% Bernoulli RV
true_p = 0.4;
pd = makedist('Binomial', 'p', true_p);
% The CLT is true for any distribution (with finite variance). Try other ones!
% pd = makedist('Exponential', 'mu', 10);
% pd = makedist('Uniform');
% sampling distribution
% generate samples and estimates
ns = [5 20 100]; % size of samples
S = 10000; % number of samples
estimates = zeros(length(ns), S);
sample_stdevs = zeros(length(ns), S);
for i = 1:length(ns)
    n = ns(i);
    for j = 1:S
        % sample from the normal RV
        x = random(pd, 1, n);
        % compute the estimate for the sample
        estimates(i, j) = mean(x);
        sample_stdevs(i, j) = std(x);
    end
end

% plot the sampling distribution
clf;
for i = 1:length(ns)
    x = estimates(i, :);
    histogram(x, 100, 'Normalization', "pdf", 'EdgeColor','none'); % if the PDF doesn't look good, increase the number of bins
    hold on
end
stem(true_p, 40, '--k', 'Linewidth', 2);
hold on
title('Sampling distribution for the mean of a Bernoulli RV') %, 'Interpreter', 'latex')
legend('sample size n = 5', 'sample size n = 20', 'sample size n = 100', 'true value', 'FontSize', 10.0);
xlabel('Estimate $\bar{X}$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);

% compute empirical Confidence Intervals
clf
alpha = 0.05;
eCIs = zeros(length(ns), 2);
for i = 1:length(ns)
    eCIs(i, 1) = quantile(estimates(i, :), alpha / 2);
    eCIs(i, 2) = quantile(estimates(i, :), 1 - alpha / 2);
    % sizes
    eCIs(i, 2) - eCIs(i, 1)
    x = estimates(i, :);
    histogram(x, 100, 'Normalization', "pdf", 'EdgeColor','none');
    hold on
end

stem(true_p, 40, '--k', 'Linewidth', 2);
hold on
stem(eCIs(1, 1), 40, '--b', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(1, 2), 40, '--b', 'Linewidth', 3,'HandleVisibility','off');
hold on
stem(eCIs(2, 1), 40, '--r', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(2, 2), 40, '--r', 'Linewidth', 3,'HandleVisibility','off');
hold on
stem(eCIs(3, 1), 40, '--y', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(3, 2), 40, '--y', 'Linewidth', 3,'HandleVisibility','off');

% compute the z-CI (Lecture 3, Eq. 12) for each sample
zCIs1 = zeros(length(ns), S); % lower edges
zCIs2 = zeros(length(ns), S); % upper edges
for i = 1:length(ns)
    zCIs1(i, :) = estimates(i, :) - norminv(1 - alpha / 2) * (estimates(i, :) .* (1 - estimates(i, :)) / ns(i)) .^ 0.5;
    zCIs2(i, :) = estimates(i, :) - norminv(alpha / 2) * (estimates(i, :) .* (1 - estimates(i, :)) / ns(i)) .^ 0.5;
    % average sizes
    mean(zCIs2(i, :) - zCIs1(i, :), 2)
end

% fraction of intervals that contain the true mean, for the various sizes n
mean(zCIs1 < true_p & true_p < zCIs2, 2) % mean over the second dimension (axis)

tCIs1 = zeros(length(ns), S); % lower edges
tCIs2 = zeros(length(ns), S); % upper edges
for i = 1:length(ns)
    tCIs1(i, :) = estimates(i, :) - tinv(1 - alpha / 2, ns(i) - 1) * sample_stdevs(i, :) / (ns(i))^0.5;
    tCIs2(i, :) = estimates(i, :) - tinv(alpha / 2, ns(i) - 1) * sample_stdevs(i, :) / (ns(i))^0.5;
    % average sizes
    mean(tCIs2(i, :) - tCIs1(i, :), 2)
end

% fraction of intervals that contain the true mean, for the various sizes n
mean(tCIs1 < true_p & true_p < tCIs2, 2) % mean over the second dimension (axis)

%% 3. Exercises
% Fruit Flies

n = 5903;
alpha = 0.01;
% point estimate
p = 637 / n

% standard error
se = (p * (1 - p) / n)^0.5

% z-CI
zCI = [p - norminv(1 - alpha / 2) * se p - norminv(alpha / 2) * se]

% Phone calls
% using a t-CI, assuming the call time distribution is normal or not too skewed:
n = 58;
alpha = 0.1;
est = 99;
std = 2576^0.5;
tCIlower = est - tinv(1 - alpha / 2, ns(i) - 1) * std / n^0.5
tCIupper = est - tinv(alpha / 2, ns(i) - 1) * std / n^0.5

% Exit Polls
n = 100;
alpha = 0.01;
% point estimate
p = 0.45
% standard error
se = (p * (1 - p) / n)^0.5

% z-CI
zCI = [p - norminv(1 - alpha / 2) * se p - norminv(alpha / 2) * se]

alpha = 2 * normcdf((p - 0.5) / se)

% z-CI
zCI = [p - norminv(1 - alpha / 2) * se p - norminv(alpha / 2) * se]

% confidence that yes will not win
1 - alpha + alpha/2

%% 4. When the t-CI fails: the mean of a Pareto distribution
% Pareto RV
theta = 1.5;
% x = rand(1, n) .^ (-1 / theta);
true_m = theta / (theta - 1)

% sampling distribution
% generate samples and estimates
ns = [5 20 100]; % size of samples
S = 1000; % number of samples
estimates = zeros(length(ns), S);
sample_stdevs = zeros(length(ns), S);
for i = 1:length(ns)
    n = ns(i);
    for j = 1:S
        % sample from the normal RV
        x = rand(1, n) .^ (-1 / theta);
        % compute the estimate for the sample
        estimates(i, j) = mean(x);
        sample_stdevs(i, j) = std(x);
    end
end

% plot the sampling distribution
clf;
for i = 1:length(ns)
    x = estimates(i, :);
    histogram(x, 1:true_m*4/100:true_m*4, 'Normalization', "pdf", 'EdgeColor','none');
    hold on
end
stem(true_m, 1, '--k', 'Linewidth', 2);
hold on
% if the PDF doesn't look good, increase the number of b
title('Sampling distribution for the mean of a Pareto RV') %, 'Interpreter', 'latex')
legend('sample size n = 5', 'sample size n = 20', 'sample size n = 100', 'true value', 'FontSize', 10.0);
xlabel('Estimate $\bar{X}$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);
% set(gca, 'XScale', 'log');
xlim([1 true_m*4]);

% compute empirical Confidence Intervals
clf
alpha = 0.05;
eCIs = zeros(length(ns), 2);
for i = 1:length(ns)
    eCIs(i, 1) = quantile(estimates(i, :), alpha / 2);
    eCIs(i, 2) = quantile(estimates(i, :), 1 - alpha / 2);
    % sizes
    eCIs(i, 2) - eCIs(i, 1)
    x = estimates(i, :);
    histogram(x, 1:true_m*4/100:true_m*4, 'Normalization', "pdf", 'EdgeColor','none');
    hold on
end

stem(true_m, 1, '--k', 'Linewidth', 2);
hold on
stem(eCIs(1, 1), '--b', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(1, 2), '--b', 'Linewidth', 3,'HandleVisibility','off');
hold on
stem(eCIs(2, 1), '--r', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(2, 2), '--r', 'Linewidth', 3,'HandleVisibility','off');
hold on
stem(eCIs(3, 1), '--y', 'Linewidth', 3,'HandleVisibility','off');
stem(eCIs(3, 2), '--y', 'Linewidth', 3,'HandleVisibility','off');
xlim([1 true_m*4]);

% compute the t-CI (Lecture 3, Eq. 13) for each sample
tCIs1 = zeros(length(ns), S);
% lower edges
tCIs2 = zeros(length(ns), S); % upper edges
for i = 1:length(ns)
    tCIs1(i, :) = estimates(i, :) - tinv(1 - alpha / 2, ns(i) - 1) * sample_stdevs(i, :) / (ns(i))^0.5;
    tCIs2(i, :) = estimates(i, :) - tinv(alpha / 2, ns(i) - 1) * sample_stdevs(i, :) / (ns(i))^0.5;
    % average sizes
    mean(tCIs2(i, :) - tCIs1(i, :), 2)
end

% fraction of intervals that contain the true value, for the various sizes n
mean(tCIs1 < true_m & true_m < tCIs2, 2) % mean over the second dimension (axis)

%% 5. Bootstrapping
% bootstrap-CI for the mean of a Normal RV

% Normal RV
mu = 0;
sigma = 1;
true_m = mu;
% observed sample
n = 5;
d = normrnd(mu, sigma, 1, n);
S = 1000;
estimates = zeros(1, S);
bootstrap_ests = zeros(1, S);
for i = 1 : S
    % sampling distribution
    x = normrnd(mu, sigma, 1, n);
    estimates(i) = mean(x);
    % this is a new sample from the Normal distrib
    % bootstrap distribution
    x = datasample(d, n); % this is sampled with replacement from the original data
    bootstrap_ests(i) = mean(x);
end
% compute the boot-CI

% plot
clf
histogram(estimates, 'normalization', 'pdf')
hold on
histogram(bootstrap_ests, 'normalization', 'pdf')
hold on
stem(true_m, 2, '--k', 'Linewidth', 3);
hold on
stem(bCI(1), 2, '-r', 'Linewidth', 4);
hold on
stem(bCI(2), 2, '-r', 'Linewidth', 4);
title('Sampling and Bootstrap distributions for the mean of a Normal RV') %, 'Interpreter', 'latex')
legend('Sampling distr', 'Bootstrap distr', 'true value', 'boot-CI', 'FontSize', 10.0);
xlabel('Estimate $\bar{X}$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);
alpha = 0.5;
bCI = quantile(bootstrap_ests, [alpha/2 1-alpha/2])

alpha = 0.05;
R = 1000;
bCIs = zeros(R, 2);
for j = 1:R
    % observed sample
    n = 5;
    d = normrnd(mu, sigma, 1, n);
    S = 1000;
    bootstrap_ests = zeros(1, S);
    for i = 1 : S
        % bootstrap distribution
        x = datasample(d, n); % this is sampled with replacement from the original data
        bootstrap_ests(i) = mean(x);
    end
    % compute the boot-CI
    bCIs(j, :) = quantile(bootstrap_ests, [alpha/2 1-alpha/2]);
end
% average size of b-CIs
mean(bCIs(:, 2) - bCIs(:, 1))

% fraction of intervals that contain the true mean, for the various sizes n
mean(bCIs(:, 1) < true_m & true_m < bCIs(:, 2), 1) % mean over the second dimension (axis)

% bootstrap-CI for the mean of a Pareto RV
% Pareto RV
theta = 1.5;
true_m = theta / (theta - 1);
% observed sample
n = 5;
d = rand(1, n) .^ (-1 / theta);
S = 1000;
estimates = zeros(1, S);
bootstrap_ests = zeros(1, S);
for i = 1 : S
    % sampling distribution
    x = rand(1, n) .^ (-1 / theta);
    estimates(i) = mean(x);
    % this is a new sample from the Pareto distrib
    % bootstrap distribution
    x = datasample(d, n); % this is sampled with replacement from the original data
    bootstrap_ests(i) = mean(x);
end
% compute the boot-CI
alpha = 0.5;
bCI = quantile(bootstrap_ests, [alpha/2 1-alpha/2])

% plot
clf
histogram(estimates, 1:true_m*3/100:true_m*3, 'normalization', 'pdf')
hold on
histogram(bootstrap_ests, 1:true_m*3/100:true_m*3, 'normalization', 'pdf')
hold on
stem(true_m, 2, '--k', 'Linewidth', 3);
hold on
stem(bCI(1), 2, '-r', 'Linewidth', 4);
hold on
stem(bCI(2), 2, '-r', 'Linewidth', 4);
title('Sampling and Bootstrap distributions for the mean of a Pareto RV') %, 'Interpreter', 'latex')
legend('Sampling distr', 'Bootstrap distr', 'true value', 'boot-CI', 'FontSize', 10.0);
xlabel('Estimate $\bar{X}$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);
xlim([1 true_m*3])

alpha = 0.05;
parameter = @(y) mean(y);
R = 1000;
bCIs = zeros(R, 2);
for j = 1:R
    % observed sample
    n = 5;
    d = rand(1, n) .^ (-1 / theta);
    S = 1000;
    ci = bootci(S, {parameter, d}, 'alpha', alpha, 'type', 'percentile');
    bCIs(j, :) = ci;
end
% average size of b-CIs
mean(bCIs(:, 2) - bCIs(:, 1))

% fraction of intervals that contain the true mean, for the various sizes n
mean(bCIs(:, 1) < true_m & true_m < bCIs(:, 2), 1) % mean over the second dimension (axis)

% bootstrap-CI for the standard deviation of a population
% observed sample
d = [176, 165, 189, 180, 172, 169, 162, 161, 183, 170];
n = length(d);
s = std(d)

S = 10000;
estimates = zeros(1, S);
bootstrap_ests = zeros(1, S);
for i = 1 : S
    % sampling distribution
    x = normrnd(mu, sigma, 1, n);
    
    estimates(i) = mean(x);
    % this is a new sample from the Normal distrib
    % bootstrap distribution
    x = datasample(d, n); % this is sampled with replacement from the original data
    bootstrap_ests(i) = std(x);
end

% compute the boot-CI
alpha = 0.05;
bCI = quantile(bootstrap_ests, [alpha/2 1-alpha/2])

% plot
clf
% histogram(estimates, 'normalization', 'pdf')
% hold on
histogram(bootstrap_ests, 'normalization', 'pdf')
hold on
stem(s, 4, '--k', 'Linewidth', 3);
hold on
stem(bCI(1), 4, '-r', 'Linewidth', 4);
hold on
stem(bCI(2), 4, '-r', 'Linewidth', 4);
title('Bootstrap distributions for the mean of a Normal RV') %, 'Interpreter', 'latex')
legend('Bootstrap distr', 'true value', 'boot-CI', 'FontSize', 16.0);
xlabel('Estimate $\bar{X} = 9.24$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);
title('Bootstrap distribution for the stdev of a Normal RV') %, 'Interpreter', 'latex')
legend('Bootstrap distr.', '$\hat{s}^* = 9.24$', '95\% boot-CIs', 'FontSize', 16.0, 'Interpreter', 'latex');
xlabel('Estimate $s$', 'Interpreter', 'latex');
ylabel('PDF');
set(gca, 'FontSize', 16.0);