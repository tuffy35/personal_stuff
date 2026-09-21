data = readmatrix('idc_lottery_test_results.csv');

sim_number = data(:,1);
cum_percent = data(:,2);

figure(1)
plot(sim_number,cum_percent)
ylabel('Cumulative Percentage of Winning')
xlabel('Simulation Number')
title('IDC Lottery Simulation when 75 IDCs Are Chosen')