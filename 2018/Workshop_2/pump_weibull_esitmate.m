pumps = xlsread('HYDPUMP.XLS');
num_pumps = 20;
pump_breakage_data=[];
for i=1:length(pumps)
pump_data = pumps(i,:);
num_weeks = pump_data(1);
num_pumps_surviving = pump_data(2);
num_pumps_broken = num_pumps-num_pumps_surviving;
num_pumps = num_pumps_surviving;
pump_breakage_data = [pump_breakage_data; num_weeks*ones(num_pumps_broken,1)];
end
[parmhat,parmci] = wblfit(pump_breakage_data, 0.05)