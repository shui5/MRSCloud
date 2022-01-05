%Run MEGAsLASER simulations
clc
clear
close all

tic
addpath /Users/steve/Documents/My_Studies/Hercules_2_Study/Code/
metabolites = {'GABA'}; % Change the targeting metabolite
freq_ppms={[1.9 7.46]}; % Change the targeting PPM
[outA,outB] = sim_MEGAsLASER_SH(metabolites,freq_ppms);
toc

% plot GABA
x_lim = [1 5];
y_lim = [-0.1 0.2];
figure(1), plot(outA.ppm,outB.specs-outA.specs,'r'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('GABA On')


