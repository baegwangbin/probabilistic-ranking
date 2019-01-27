% Answer for Question B
% Run Question_B_1.m before running this

seed_1 = 1;
seed_2 = 2;
seed_3 = 3;

figure
ax1 = subplot(3,1,1);
fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
hold on; plot(1:200, P_B(:,seed_1), 'r-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,200])
ylim([0,3])
ylabel('seed = 1', 'FontSize', 15,'FontWeight','bold');   
title('Changing Random Seed for Same Player', 'FontSize', 20,'FontWeight','bold')
   
ax2 = subplot(3,1,2);
fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
hold on; plot(1:200, P_B(:,seed_2), 'g-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,200])
ylim([0,3])
ylabel('seed = 2', 'FontSize', 15,'FontWeight','bold');   

ax3 = subplot(3,1,3);
fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
hold on; plot(1:200, P_B(:,seed_3), 'b-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,200])
ylim([0,3])
xlabel('Gibbs Iteration', 'FontSize', 15,'FontWeight','bold');
ylabel('seed = 3', 'FontSize', 15,'FontWeight','bold');   
