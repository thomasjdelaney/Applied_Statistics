%% Workshop 10
%% Question 1.1
% reading data
delimiterIn = ',';
A = importdata('3dprinters.txt', delimiterIn);

% plotting data
clf
hist(A.data(:,1))
xlabel('failure (yes/no)')
ylabel('frequency')

% we will look at unique temperatures recorded.
tmp = unique(A.data(:,2));
fraction = zeros(0,length(tmp));
for i = 1:length(tmp)
    idx = find(A.data(:,2)==tmp(i));
    fraction(i) = sum(A.data(idx,1))/length(idx);
end
clf
plot(tmp,fraction,'.')
xlabel('temperature (degrees Celsius)')
ylabel('proportion of successful prints')

%% Question 1.2
% we start by defining our variables:
response = A.data(:,1);
temperature = A.data(:,2);
% need 2 dummy variables for polymer type. Level 0 is absorbed into the intercept.
dummy1 = zeros(length(response),1);
dummy2 = zeros(length(response),1);
dummy1(find(A.data(:,3)==1)) = 1;
dummy2(find(A.data(:,3)==2)) = 1;
 
% Now define values for parameters that we use at the start of the optimisation.
% We have 4 parameters:
x0 = [0 0 0 0];

% check that computing the negative log-likelihood
% usning 'negloglikGLM.m' works:
negloglikGLM([1 2 0 1],1,2,0,0) % uses beta=[1 2 0 1].
% response=1, temperature=20, dummy1=0, dummy=1
 
% maximum likelihood fit:
fun = @(s)negloglikGLM(s,response',temperature',...
dummy1',dummy2');
options = optimset('MaxFunEval',10000000,'MaxIter',10000);
s = fminsearch(fun,x0,options); % here we performe the optimisation
% (i.e. find parameters for which negloglik is minimal)
 
% The fitted parameters are stored in 's':
disp(s)

%% Question 1.3
% create a table for the data:
response = nominal(response);
polymer = A.data(:,3);
polymer = nominal(polymer);
data = table(response,temperature,polymer,'VariableNames',{'success','temperature','polymer'});
m1 = fitglm(data,'success~temperature+polymer','Distribution','binomial');

% setting up the service duration interval
% we want to plot model fit lines for:
w = linspace(min(temperature),max(temperature));
% plot the data and fit lines:
clf
plot(tmp,fraction,'.')
xlabel('temperature (degrees Celsius)')
ylabel('proportion of successful prints')
line(w,feval(m1,w,'0'),'Color','b','LineWidth',1)
line(w,feval(m1,w,'1'),'Color','g','LineWidth',1)
line(w,feval(m1,w,'2'),'Color','r','LineWidth',1)

%% Question 1.4
% residual plots
% Q-Q plot to check normality
subplot(2,2,1)
plotResiduals(m1,'probability','ResidualType','Deviance')% residuals versus fitted values
subplot(2,2,2)
plotResiduals(m1,'fitted','ResidualType','Deviance')% auto-correlation (via lagged residuals)
subplot(2,2,3)
plotResiduals(m1,'lagged','ResidualType','Deviance')% leverage plot
subplot(2,2,4)
plotDiagnostics(m1)

% fit a constant model to the data:
m0 = fitglm(data,'success~1','Distribution','binomial');
[h pvalue stat] = lratiotest(m1.LogLikelihood,m0.LogLikelihood,3);

%% Question 2.1
% read data:
delimiterIn = ' ';
A = importdata('ship_damage.txt',delimiterIn);
service = A.data(:,1);
damages = A.data(:,2); 
type = A.textdata(2:35,2);
type = nominal(type);
 
clf
subplot(1,2,1)
hist(damages)
title('distribution of damage incident counts')
xlabel('damage incident count')
ylabel('frequency')
subplot(1,2,2)
gscatter(A.data(:,1),A.data(:,2),A.textdata(2:35,2),'bgr','x.o*')
xlabel('months in service'); ylabel('number of damage incidents');
title('damage in different ship types')

%% Question 2.2
data = table(damages,service,type,'VariableNames',{'damage_count','service_duration','ship_type'});
m1 = fitglm(data,'damage_count~service_duration+ship_type','Distribution','Poisson');

% setting up the service duration interval
% we want to plot model fit lines for:
w = linspace(min(service),max(service));
 
% plot the data and fit lines:
clf
gscatter(service,damages,A.textdata(2:35,2),'bgryc','x.o*');
xlabel('months in service'); ylabel('number of damage incidents');
title('damage in different ship types');
line(w,feval(m1,w,'A'),'Color','b','LineWidth',1);
line(w,feval(m1,w,'B'),'Color','g','LineWidth',2);
line(w,feval(m1,w,'C'),'Color','r','LineWidth',1);
line(w,feval(m1,w,'D'),'Color','y','LineWidth',1);
line(w,feval(m1,w,'E'),'Color','c','LineWidth',1);

% residual plots
% Q-Q plot to check normality
subplot(2,2,1)
plotResiduals(m1,'probability','ResidualType','Deviance')% residuals versus fitted values
subplot(2,2,2)
plotResiduals(m1,'fitted','ResidualType','Deviance')% auto-correlation (via lagged residuals)
subplot(2,2,3)
plotResiduals(m1,'lagged','ResidualType','Deviance')% leverage plot
subplot(2,2,4)
plotDiagnostics(m1)