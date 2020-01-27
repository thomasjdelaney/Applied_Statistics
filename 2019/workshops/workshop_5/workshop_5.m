%% Workshop 5
%% Question 1 (a)
%%
% 
x = [15.1 14.6 12.0 19.2 16.1 15.5 11.3 18.7 17.1 17.2]; % input vector
parameter = @(y)mean(y);
ci=bootci(3000,{parameter, x}, 'alpha', 0.05,'type', 'percentile');
[ci,bootstat] = bootci(3000,{parameter, x}, 'alpha', 0.05, 'type', 'percentile');
histogram(bootstat)
bootstrap_se = sqrt(sum((bootstat - mean(x)).^2)/(3000-1));
percentiles = prctile(bootstat, [5, 95]);
% similar to ci above, but percentiles CI is tighter

[tci,tbootstat] = bootci(3000,{parameter, x}, 'alpha', 0.05, 'type', 'student'); % widest
[nci,nbootstat] = bootci(3000,{parameter, x}, 'alpha', 0.05, 'type', 'normal'); % wider

%% Question 1 (b)
%%
samples = lognrnd(2, 0.5, 10, 3000);
samples_means = mean(samples);
random_index = randi([1 3000]);
random_sample = samples(:, random_index);

[pci, pbootstat]=bootci(3000,{parameter, random_sample}, 'alpha', 0.05,'type', 'percentile');
[sci, sbootstat]=bootci(3000,{parameter, random_sample}, 'alpha', 0.05,'type', 'student');
[nci, nbootstat]=bootci(3000,{parameter, random_sample}, 'alpha', 0.05,'type', 'normal');

