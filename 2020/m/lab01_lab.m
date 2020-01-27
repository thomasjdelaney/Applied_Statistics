%% lab01_lab.m
%% Random Variables, probability distributions, and sampling
%% 1. Random Variables
% Binomial distribution
N = 10;
p = 0.5;

pd = makedist('Binomial', 'N', N, 'p', p);

% Theoretical PMF
outcomes = 0:1:10;
binom_pmf = pdf(pd, outcomes);
clf;
plot(outcomes, binom_pmf, 'o-b', 'Linewidth', 2);
hold on;
binom_mean = mean(pd);
stem(binom_mean, max(binom_pmf), '--r', 'Linewidth', 3);
legend('theoretical PDF', 'mean', 'Location', 'northeast');
title('PDF - Binomial');
xlabel('outcomes');
ylabel('PDF');

% Generate a random sample
rng('default');
S = 10; % sample size
small_sample = random(pd, 1, S);
rng('default');
binornd(N, p, 1, S);

% Empirical PMF
binom_epmf = histogram(small_sample, 'normalization', 'pdf');
binom_epmf.Values;

large_sample = random(pd, 1, 10000);
histogram(large_sample, 'BinEdges', (0.5:1:N + 0.5), 'normalization', 'pdf');

% Theoretical CDF
binom_cdf = cdf(pd, outcomes);
clf; 
stairs(outcomes, binom_cdf, 'o-r', 'LineWidth', 2);
hold on;
legend('theoretical CDF', 'Location', 'southeast');
title('CDF - Binomial');
xlabel('outcomes');
ylabel('CDF');

% Empirical CDF
[binom_ecdf, x] = ecdf(small_sample);
stairs(x, binom_ecdf, 'o-');
legend('theoretical CDF','empirical CDF');

[binom_ecdf, x] = ecdf(large_sample);
stairs(x, binom_ecdf, 's--k', 'LineWidth', 2)
legend('theoretical CDF','small sample', 'large sample')

% Quantile function
icdf(pd, 0.1);
icdf(pd, 0.2);
% evaluate at more points within 0 and 1
pp = 0: 0.01: 1;icdf(pd, 0.1);
icdf(pd, 0.2);
% evaluate at more points within 0 and 1
pp = 0: 0.01: 1;

clf;
stairs(pp, icdf(pd, pp), '-');
hold on;
legend('theoretical quantile function', 'Location', 'southeast');
title('Inverse CDF - Binomial')
xlabel('outcomes')
ylabel('CDF')

quantile(large_sample, 0.1)
quantile(large_sample, 0.2)
% evaluate at more points within 0 and 1
stairs(pp, quantile(large_sample, pp), '--')
legend('theoretical quantile function', 'large sample empirical quantile function', 'Location','southeast')

% Q-Q plot
clf;

theor_q = icdf(pd, pp);
empir_q = quantile(small_sample, pp);
empir_q_ls = quantile(large_sample, pp);

plot(theor_q, empir_q, 'o');
hold on;
plot(theor_q, empir_q_ls, 's', 'MarkerSize', 10);
hold on;
plot(theor_q, theor_q, '-', 'HandleVisibility', 'off');
legend('small sample', 'large sample', 'Location', 'southeast');
title('Q-Q plot - Binomial') % title for plot
xlabel('theoretical quantiles') % x-axis label
ylabel('empirical quantiles') % y-axis label

% Common distributions
% 1,2.  Uniform, continuous, [a,b], (a,b)
%       Exponential, continuous, [0, infinity), rate
%       Multinomial, discrete, integer array summing to n, p_i summing to 1
%       Gamma, continuous, (0, infinity), (shape, rate)
%       Normal, continuous, (-infinity, infinity), (mean, variance)

%% 2. Estimate the probability of events

% 1. Exponential
pd = makedist('Exponential', 'mu', 4);
% cdf(pd, 6) - cdf(pd, 5);

% 2. Normal
sigma = 1; % try different values of sigma, the answer will always be the same.
pd = makedist('Normal', 'mu', 0, 'sigma', sigma);
% cdf(pd, sigma) - cdf(pd, -sigma);

% 3. Uniform
pd = makedist('Uniform', 'lower', 0, 'upper', 60);
% cdf(pd, 20) - cdf(pd, 5);
% cdf(pd, 45) - cdf(pd, 30);
% cdf(pd, 70) - cdf(pd, 55)

%% 3. Compute the empirical CDF from data
data = [0.9838 1.1704 6.8989 1.7893 1.5168 0.1136 6.3209 1.6271 3.8187 3.5633];
% sort data in ascending order, x
x = sort(data);
% create array y of increasing integers of size length(x)
y = (1: 1: length(x));
% divide y by length(x)
y = y ./ length(x);
% Matlab's CDF
[y_ecdf, x_ecdf] = ecdf(data);

clf;
stairs(x, y, '-k', 'LineWidth', 2);
hold on
stairs(x_ecdf, y_ecdf, '--r', 'LineWidth', 2)
legend('our ECDF', 'Matlab ECDF', 'Location', 'southeast')
title('Empirical CDF') % title for plot
xlabel('outcomes') % x-axis label
ylabel('ECDF') % y-axis label

%% 4. Compute the empirical PDF from data
data = [0.9838 1.1704 6.8989 1.7893 1.5168 0.1136 6.3209 1.6271 3.8187 3.5633];
bin_edges = 0: 3: 12;
% compute bins' mid points
x = (bin_edges(1: end-1) + bin_edges(2: end)) / 2;
% count the number of data points that fall in each bin
[c, edges] = histcounts(data, bin_edges);
% normalise
y = zeros(1, length(c));
for i = 1:length(c)
    y(i) = c(i) / (length(data) * (bin_edges(i+1) - bin_edges(i)));
end

clf;
histogram(data, 'BinEdges', bin_edges, 'normalization', 'pdf');
hold on
bar(x, y)
legend('Matlab PDF', 'our PDF', 'Location', 'northeast')
title('PDF') % title for plot
xlabel('outcomes') % x-axis label
ylabel('PDF') % y-axis label

%% 5. Inverse Probability integral transform sampling
% generate uniform random numbers
unif = makedist('Uniform', 'lower', 0, 'upper', 1);
u_samples = random(unif, 1000);

% compute the inverse CDF of the desired distribution
pd = makedist('Normal', 'mu', 5, 'sigma', 3);
x_samples = icdf(pd, u_samples);

% PDFs
clf;
histogram(x_samples, 'normalization', 'pdf');
hold on;
x = -10: 0.1: 20;
th_pdf_x = pdf(pd, x);
plot(x, th_pdf_x, '-r', 'LineWidth', 3)
legend('sample', 'theory', 'Location', 'northwest')
title('PDF Normal') % title for plot
xlabel('x') % x-axis label
ylabel('PDF') % y-axis label

% sum of two dices
outcomes = (2 : 1: 12);
th_pdf = [1 2 3 4 5 6 5 4 3 2 1] / 36;
th_cdf = cumsum(th_pdf);

% generate uniform random numbers
unif = makedist('Uniform', 'lower', 0, 'upper', 1);
u_samples = random(unif, 1000);

x_samples = zeros(1, length(u_samples));
for i = 1: length(u_samples)
    u = u_samples(i);
    % Inverse PIT:
    x_samples(i) = outcomes(min(find(th_cdf > u)));
end

% PDFs
clf; % clear the figure
histogram(x_samples, 'normalization', 'pdf');
hold on
% Theoretical PDF
plot(outcomes, th_pdf, '-r', 'LineWidth', 2)
legend('sample', 'theory', 'Location', 'northwest')
title('PDF Normal') % title for plot
xlabel('x') % x-axis label
ylabel('PDF') % y-axis label

%% 6. Functions of random variables and sampling
% Sample
unif = makedist('Uniform', 'lower', 0, 'upper', 1);
u_samples = random(unif, 1000);
x_samples = u_samples .^ (2/3);

% PDFs
clf;
histogram(x_samples, 'normalization', 'pdf');
hold on;
x = 0: 0.001: 1;
th_pdf_x = 3/2 .* x .^ (1/2);
plot(x, th_pdf_x, '-r', 'LineWidth', 3)
legend('sample', 'theory', 'Location', 'northwest')
title('PDF X') % title for plot
xlabel('x') % x-axis label
ylabel('PDF') % y-axis label

% Transformation Y = e^X
y_samples = exp(x_samples);
clf;
histogram(y_samples, 'normalization', 'pdf');
hold on;
y = exp(0): 0.001: exp(1); % remeber to transform the domain as well
th_pdf_y = 3/2 .* log(y) .^ (1/2) ./ abs(y);
plot(y, th_pdf_y, '-r', 'LineWidth', 3)
legend('sample', 'theory', 'Location', 'northwest')
title('PDF Y') % title for plot
xlabel('y') % x-axis label
ylabel('PDF') % y-axis label

% Transformation Z = e^{-X}
y_samples = exp(- x_samples);
clf;
histogram(y_samples, 'normalization', 'pdf');
hold on;
z = exp(-1): 0.001: exp(-0); % remeber to transform the domain as well
th_pdf_z = 3/2 .* (-log(z)) .^ (1/2) ./ abs(z);
plot(z, th_pdf_z, '-r', 'LineWidth', 3)
legend('sample', 'theory', 'Location', 'northeast')
title('PDF Z') % title for plot
xlabel('z') % x-axis label
ylabel('PDF') % y-axis label

%% 7. Donations
% Distribution of one donation
pd = makedist('Exponential', 'mu', 5);

% Samples of N donations
S = 1000;
N = 10;
y_sample = sum(random(pd, N, S), 1);

% PDF
clf;
histogram(y_sample, 'normalization', 'pdf');
title('PDF Y') % title for plot
xlabel('y') % x-axis label
ylabel('PDF') % y-axis label

% Empirical Mean
mean(y_sample);

% Probability to collect more than 60 pounds
T = 60;
Ss = round(logspace(2, 7, 40));
% Ss = 1000:1000:300000;
CCDFyT = zeros(1, length(Ss));
for i = 1:length(Ss)
    S = Ss(i);
    y_sample = sum(random(pd, N, S), 1);
    CCDFyT(i) = sum(y_sample > T) / S;
end
clf;
plot(Ss, CCDFyT, '-o')
title('$CCDF_Y(T)$ vs $S$', 'Interpreter', 'latex') % title for plot
xlabel('S, sample size') % x-axis label
ylabel('$CCDF_Y(T)$', 'Interpreter', 'latex') % y-axis label
set(gca, 'XScale', 'log')