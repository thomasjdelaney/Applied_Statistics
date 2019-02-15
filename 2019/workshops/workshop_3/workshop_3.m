%% Applied statistics workshop 3
%% Question 1 (a)
m = 2500;
x = [2610, 2750, 2420, 2510, 2540, 2490, 2680];
num_samples = length(x);
mean_x = mean(x);
std_x = std(x);
alpha = 0.1;
tail = 'right'; % sample mean should be greater than or equal to m
[h, p, ci, stats] = ttest(x, m, alpha, tail); % [1, 0.0759, 2.5088e^3, stats]
tv = tinv(0.1, 6); % rejection region of t-dist w/ 6 degrees of freedom
power = sampsizepwr('t', [m std_x], mean_x, [], 7, 'Tail', 'right', 'Alpha', alpha);
sample_size = sampsizepwr('t', [m std_x], mean_x, 0.9, [], 'Tail', 'right', 'Alpha', alpha);
nn = 1:70;
plot_power = sampsizepwr('t', [m std_x], mean_x, [], nn, 'Alpha', 0.1);
plot_power_alpha = sampsizepwr('t', [m std_x], mean_x, [], nn, 'Alpha', 0.05);
plot(plot_power);
hold on
plot(plot_power_alpha, 'r*');
title('Power vs Sample Size');
xlabel('Sample Size');
ylabel('Power');

%% Question 1(b)
n = 40;
m = 48;
num_samples = 40;
mean_x = 60;
std_x = 41;
alpha = 0.01;
tail = 'right';
tv = (m - mean_x)/(std_x/sqrt(num_samples));
tstat = tinv(0.1, 40);

[h, p, ci, z_val] = ztest(mean_x, m, std_x, 'Alpha', 0.1, 'Tail', 'right');
z_value = (60 - 48)/(41/sqrt(n));
z_stat = norminv(0.1, 48, 41);
%% compare tv with the t-value?

%% Question 1(c)
tooth_depth = xlsread('CHEEKTEETH.xls');
num_samples = length(tooth_depth);
mean_x = mean(tooth_depth);
std_x = std(tooth_depth);
m = 15;
[h, p, ci, stats] = ttest(tooth_depth, m); % [1, 0.049, [15.5, 17.5], stats]

%% Question 2(a)
num_samples = 351;
num_successes = 144;
phat = num_successes/num_samples;
p0 = 0.35;
std_dev = sqrt(p0*(1 - p0));
Z = (phat - p0)/sqrt(p0*(1 - p0)/num_samples);
rejection_region = norminv(0.025); % -0.0242
% Z < 2.36 and so H0 is rejected
[h, p, ci, stats] = ztest(phat, p0, std_dev/sqrt(num_samples-1));

%% Question 3
antimony=xlsread('ANTIMONY.XLS');
strength=antimony(:,3);
mean_s = mean(strength);
std_s = std(strength);
X=min(strength):0.1:max(strength);
Y=normcdf(X, mean_s, std_s);
cdf=[X' Y'];
H = kstest(strength, cdf);
normplot(strength);
hist(strength);
