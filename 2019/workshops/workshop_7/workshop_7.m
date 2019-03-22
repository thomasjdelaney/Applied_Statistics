%% Workshop 7
%% Preliminary analysis

delimiterIn = ',';
A = importdata('APPLIANCES2.csv',delimiterIn);
 
%% scatterplots of data
appliances = A.data(:,1);
temperature = A.data(:,2);
pressure = A.data(:,3);
humidity = A.data(:,4);
windspeed = A.data(:,5);
 
clf
subplot(2,2,1)
plot(temperature,appliances,'.')
xlabel('temperature (Celsius)')
ylabel('appliance energy use (Wh)')
 
subplot(2,2,2)
plot(pressure,appliances,'.')
xlabel('pressure (mm hg)')
ylabel('appliance energy use (Wh)')
 
subplot(2,2,3)
plot(humidity,appliances,'.')
xlabel('humidity (Percent)')
ylabel('appliance energy use (Wh)')
 
subplot(2,2,4)
plot(windspeed,appliances,'.')
xlabel('windspeed (m/s)')
ylabel('appliance energy use (Wh)')

% create a table for the data:
data = table(temperature,pressure,humidity,windspeed,appliances,'VariableNames',{'temperature','pressure','humidity','windspeed','appliances'});
% show first few rows of table
data(1:5,:)

% start with a full model that includes all predictors:
m1 = fitlm(data, 'appliances~temperature+pressure+humidity+windspeed')