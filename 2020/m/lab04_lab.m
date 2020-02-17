%% lab04_lab.m
%% Lab 4: Hypothesis testing
%% 1. z-test for a proportion: Exit polls revisted

% true value
q0 = 0.5;
mu = q0;
n = 100;
sigma = (q0*(1-q0)/n)^0.5;
observation = 0.45;
% left-tailed p-value
pval = normcdf(observation, mu, sigma)

% significance level
alpha = 0.01;
% check
[h, p, ci, zval] = ztest(observation, mu, sigma, 'alpha', alpha, 'tail', 'left')

% plot
clf;
xx = 0.3:0.001:0.7;
plot(xx, normpdf(xx, mu, sigma), '-', 'LineWidth', 2);
hold on
stem(observation, 8, '--r', 'LineWidth', 2);
hold on
stem(norminv(alpha, mu, sigma), 8, '-.k', 'LineWidth', 2);
hold on
stem(mu, 8, ':g', 'LineWidth', 2);
legend('Sampling distrib', 'Observation', 'Threshold significant obs', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Sampling distribution')

% n needed for 0.45 to be 1% significant
n = ceil(((q0*(1-q0))^0.5/(observation - mu) * norminv(alpha))^2)

sigma = (q0*(1-q0)/n)^0.5;
observation = 0.45;
% left-tailed p-value
pval = normcdf(observation, mu, sigma)

% significance level
alpha = 0.01;
% plot
clf;
xx = 0.3:0.001:0.7;
plot(xx, normpdf(xx, mu, sigma), '-', 'LineWidth', 2);
hold on
stem(observation, 20, '--r', 'LineWidth', 2);
hold on
stem(norminv(alpha, mu, sigma), 20, '-.k', 'LineWidth', 2);
hold on
stem(mu, 20, ':g', 'LineWidth', 2);
legend('Sampling distrib', 'Observation', 'Threshold significant obs', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Sampling distribution')

%% 2. Paired t-test
m_before = csvread("m_before.csv");
m_after = csvread("m_after.csv");
Dm = m_after - m_before;
observation = mean(Dm)

clf;
histogram(Dm)

n = length(Dm);
mu = 0;
sigma = std(Dm) / n^0.5;
% test statistic
tstat = (observation - mu)/ sigma

% degrees pf freedom
r = n - 1

% two-tailed p-value
pval = 2 * (1 - tcdf(tstat, r))

% significance threshold
alpha = 0.001;
% check
[h, p, ci, stats] = ttest(Dm, mu, 'alpha', alpha, 'tail', 'both')

% plot
clf;
xx = -5:0.01:5;
plot(xx, tpdf(xx, r), '-', 'LineWidth', 2);
hold on
stem(tstat, 0.5, '--r', 'LineWidth', 2);
hold on
stem(tinv(1 - alpha/2, r), 0.5, '-.k', 'LineWidth', 2);
hold on
stem(tinv(alpha/2, r), 0.5, '-.k', 'LineWidth', 2, 'HandleVisibility','off');
hold on
stem(0, 0.5, ':g', 'LineWidth', 2);
legend('Sampling distrib', 'Observation', 'Threshold significant obs', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Sampling distribution')

%% 3. A/B test
m_new = csvread("m_new.csv");
m_old = csvread("m_old.csv");
observation = mean(m_new) - mean(m_old)

mx = mean(m_new);
sx = std(m_new);
nx = length(m_new);
my = mean(m_old);
sy = std(m_old);
ny = length(m_old);
% degrees of freedom
r = nx + ny - 2;
% test statistic
tstat = (mx - my) / ((sx^2 *(nx-1) + sy^2*(ny-1))/r *(1/nx + 1/nx))^0.5

% p value
pval = 2 * (1 - tcdf(tstat, r))

% significance threshold
alpha = 0.001;
% check
[h, p, ci, stats] = ttest2(m_old, m_new, 'alpha', alpha, 'tail', 'both')

% plot
clf;
xx = -5:0.01:5;
plot(xx, tpdf(xx, r), '-', 'LineWidth', 2);
hold on
stem(tstat, 0.5, '--r', 'LineWidth', 2);
hold on
stem(tinv(1 - alpha/2, r), 0.5, '-.k', 'LineWidth', 2);
hold on
stem(tinv(alpha/2, r), 0.5, '-.k', 'LineWidth', 2, 'HandleVisibility','off');
hold on
stem(0, 0.5, ':g', 'LineWidth', 2);
legend('Sampling distrib', 'Observation', 'Threshold significant obs', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Sampling distribution')

%% 4. Pregnancy test
mu0 = 10;
std0 = 2;
mu1 = 18;
std1 = 3;
% value of hCG concentration such that the probability to measure
% a higher value in a non-pregnant female is 0.1.
alpha = 0.1;
observation = norminv(1 - alpha, mu0, std0)

% probability to commit a type II error (failure to detect pregnancy)
% if significance level is 0.1
beta = normcdf(observation, mu1, std1)

% probability to commit a type I error (erroneous detection of pregnancy)
% if we declare a pregnancy with a hCG concentration of 13 or higher
alpha = 1 - normcdf(13, mu0, std0)

% plot
clf
x = -5:0.1:35;
plot(x, normpdf(x, mu0, std0), '-', 'LineWidth', 2);
hold on
plot(x, normpdf(x, mu1, std1), '-', 'LineWidth', 2);
hold on
stem(observation, 0.25, '--r', 'LineWidth', 2);
hold on
r = observation:0.1:35;
inBetween = [normpdf(r, mu0, std0), fliplr(zeros(1, length(r)))];
fill([r fliplr(r)], inBetween, 'g');
hold on
r = 0:0.1:observation;
inBetween = [normpdf(r, mu1, std1), fliplr(zeros(1, length(r)))];
fill([r fliplr(r)], inBetween, 'y');
legend('H_0', 'H_1', 'answer 1', 'prob type I error', 'prob type II error')

%% 5. Pregnancy test 2 (confusion matrix)
