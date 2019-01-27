% Answer for Question A
% Gibbs sampling for player skills

figure
ax1 = subplot(3,1,1);
plot(1:iteration, P_1, 'r-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,1100])
ylabel('Player 16', 'FontSize', 15,'FontWeight','bold');   
title('Player Skill vs. Gibbs Iteration', 'FontSize', 20,'FontWeight','bold')
   
ax2 = subplot(3,1,2);
plot(1:iteration, P_2, 'g-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,1100])
ylabel('Player 69', 'FontSize', 15,'FontWeight','bold');   

ax3 = subplot(3,1,3);
plot(1:iteration, P_3, 'b-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,1100])
xlabel('Gibbs Iteration', 'FontSize', 15,'FontWeight','bold');
ylabel('Player 107', 'FontSize', 15,'FontWeight','bold');   
