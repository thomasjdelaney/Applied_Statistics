%% lab05_lab.m
%% Lab 5: Nonparametric and permutation methods
%% 1. Kolmogorov-Smirnov test for one sample
% 1.
pd = makedist('Exponential', 'mu', 1);

% 4.
% pd = makedist('Normal');

n = 100;
r = 1000;

% compute the K-S test statistic
Ds = zeros(1,r);
for i = 1:r
    x = random(pd, n, 1);
    [f, x_values] = ecdf(x);
    % K-S test statistic: sup of distances
    Ds(i) = max( max(abs(f - cdf(pd, x_values))), max(abs(f(1:end-1,:) - cdf(pd, x_values(2:end,:)) )) );
end

% 2.
K = Ds * n^0.5;

% plot Kolmogorov PDF
clf
histogram(K, 'normalization', 'pdf')
hold on
k = 0:0.01:3;
plot(k, kolmpdf(k), '-', 'Linewidth', 3);

% 3.
K_accurate = Ds * n^0.5 + 1/(6 * n^0.5) + (Ds * n^0.5 - 1)/4/n;
% plot Kolmogorov PDF
% clf
histogram(K_accurate, 'normalization', 'pdf')
hold on
k = 0:0.01:3;
plot(k, kolmpdf(k), '-', 'Linewidth', 3);

%% One-sample K-S test on a Exponential distribution
% 1.
n = 100;
pd = makedist('Exponential', 'mu', 1);
x = random(pd, n, 1);

% 5.
% pd2 = makedist('Exponential', 'mu', 1.2);
% x = random(pd2, n, 1);

% plot
[f, x_values] = ecdf(x);

clf
J = stairs(x_values,f);
hold on;
xx = 0:0.1:10;
K = plot(xx, cdf(pd, xx),'r--');

set(J,'LineWidth',2);
set(K,'LineWidth',2);
legend([J K],'Empirical CDF','Theoretical CDF','Location','SE');
set(gca, 'FontSize', 16.0);
xlabel('x')
ylabel('CDF(x)')

% 2.
% K-S test statistic: supremum of distances
D = max( max(abs(f - cdf(pd, x_values))), max(abs(f(1:end-1,:) - cdf(pd, x_values(2:end,:)) )) )

% 3.
% p-value
1 - kolmcdf(D * n^0.5)

% better approximation of the p-value
1 - kolmcdf(D * n^0.5 + 1/(6 * n^0.5) + (D * n^0.5 - 1)/4/n)

%4.
% check
[h, p, ksstat] = kstest(x, 'CDF', pd, 'alpha', 0.01)

%% 2. Kolmogorov-Smirnov test for two samples
% Two-sample K-S test on a Exponential distribution
% 1.
pd = makedist('Exponential', 'mu', 1);
n = 100;
x = random(pd, n, 1);
[fx, x_values] = ecdf(x);

% 2.
m = 150;
y = random(pd, m, 1);
[fy, y_values] = ecdf(y);

% plot
clf
J = stairs(x_values,fx);
hold on;
K = stairs(y_values,fy, '-r');

set(J,'LineWidth',2);
set(K,'LineWidth',2);
legend([J K],'Sample 1','Sample 2','Location','SE');
set(gca, 'FontSize', 16.0);
xlabel('x')
ylabel('CDF(x)')

% 3.
[h,p,ks2stat] = kstest2(x, y, 'alpha', 0.01)

% compare to the p-value from the formula for large samples given in the lecture
pval = 1 - kolmcdf(ks2stat * (n*m/(n+m))^0.5)

% 4.
r = 10000;
n = 100;
m = 150;

n_rejections = 0;
for i = 1:r
    x = random(pd, n, 1);
    y = random(pd, m, 1);
    [h,p,ks2stat] = kstest2(x, y, 'alpha', 0.01);
    n_rejections = n_rejections + h;
end

typeIerr = n_rejections / r

% 5.
pd2 = makedist('Exponential', 'mu', 2);
r = 10000;
n = 100;
m = 150;
n_rejections = 0;
for i = 1:r
    x = random(pd, n, 1);
    y = random(pd2, m, 1);
    [h,p,ks2stat] = kstest2(x, y, 'alpha', 0.01);
    n_rejections = n_rejections + h;
end
% the probability of type II errors is equal to the fraction of times we
% didn't reject the null hypothesis (which is false)
typeIIerr = (r - n_rejections) / r

%% 3. Permutation test for equal means of two populations (A/B test)
g2 = csvread("m_new.csv");
g1 = csvread("m_old.csv");
% observed test statistic
observation = mean(g2) - mean(g1)

% permutation test
% randomly reassigning individuals to the two groups
n1 = length(g1);
n2 = length(g2);
g12 = cat(2, g1, g2);

S = 10000;
tstat = zeros(1, S);
for i = 1:S
    r = randperm(n1 + n2);
    gr = g12(r);
    gr1 = gr(1:n1);
    gr2 = gr(n1+1:n1+n2);
    m = mean(gr2) - mean(gr1);
    tstat(i) = m;
end
% right-tailed p-value
pval = mean(tstat >= observation)

clf
histogram(tstat, 'normalization', 'pdf');
hold on
stem(observation, 0.1, '--r', 'LineWidth', 4);
hold on
stem(0.0, 0.1, ':g', 'LineWidth', 4);
legend('Permutation distrib', 'Observation', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Permutation distribution')

%% 4. Permutation test for matched pairs
g1 = transpose(csvread("m_before.csv"));
g2 = transpose(csvread("m_after.csv"));
Dm = g2 - g1;
observation = mean(Dm)

% permutation test
% randomly swap each user's visit times to the new and old versions
n1 = length(g1);
n2 = length(g2);
g12 = cat(1, g1, g2);
S = 10000;
tstat = zeros(1, S);
for i = 1:S
    r = randi(2, 1, n1) - 1;
    % rm is a random matrix where each row can be either [-1 1] or [1 -1] with equal probability
    rm = (cat(1, r, 1 - r) - 0.5) * 2;
    gr = g12 .* rm;
    m = mean(gr(2, :) + gr(1, :));
    tstat(i) = m;
end
% right-tailed p-value
pval = mean(tstat >= observation)

clf
histogram(tstat, 'normalization', 'pdf');
hold on
stem(observation, 0.1, '--r', 'LineWidth', 4);
hold on
stem(0.0, 0.1, ':g', 'LineWidth', 4);
legend('Permutation distrib', 'Observation', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Permutation distribution')

%% 5. Permutation test for a relationship
% load temperature data
oxford = readtable('oxford_data.txt', 'Delimiter', '|');
oxford = str2double(replace(table2array(oxford(:,[3,4])), '*', ''));

% midrange temperatures
temp = mean(oxford, 2);

% numer of months since January 1853
time = transpose(1:length(temp));

% plot midrange temperature time series
clf
plot(time, temp);
set(gca, 'FontSize', 16.0);
xlabel('Months since Jan 1853')
ylabel('Midrange temperature [C]')

% use Pearson's correlation coefficient as test statistic
cc = corrcoef(time, temp, 'Rows', 'complete');
observation = cc(1, 2)

% permutation test
% randomly shuffle the time variable
S = 10000;
tstat = zeros(1, S);
for i = 1:S
    time = randperm(length(temp));
    cc = corrcoef(time, temp, 'Rows','complete');
    tstat(i) = cc(1, 2);
end
% right-tailed p-value
pval = mean(tstat >= observation)

% plot
clf
histogram(tstat, 'normalization', 'pdf');
hold on
stem(observation, 20, '--r', 'LineWidth', 4);
hold on
stem(0.0, 20, ':g', 'LineWidth', 4);
legend('Permutation distrib', 'Observation', 'H_0')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Permutation distribution')

%% 6. Constrained permutation test
% let's create the data reported in the table
% vector of student genders (20 females, 30 males)
% 1 female, 0 male
gender = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% total number of students
n = length(gender)

% number of female students
n_females = sum(gender)

% vector of exam outcomes (35 pass, 10 fail)
% 1 fail, 0 pass
result = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% number of failures
n_fail = sum(result)

% number of females who failed the exam
observation = sum(gender .* result)

% permutation test
% randomly reassign the exam outcomes
S = 10000;
tstat = zeros(1, S);
for i = 1:S
% number of females who failed the exam
    tstat(i) = sum(gender(randperm(n)) .* result);
end
% two-tailed p-value
2 * mean(tstat <= observation)

% plot
clf
histogram(tstat, 'normalization', 'pdf');
hold on
stem(observation, 0.25, '--r', 'LineWidth', 4);
legend('Permutation distrib', 'Observation')
set(gca, 'FontSize', 16.0);
xlabel('Test statistic')
ylabel('Permutation distribution')