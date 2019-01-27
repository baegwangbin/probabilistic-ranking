% Answer for Question A
% Plot sampled skills (burn-in & thinning)

figure
ax1 = subplot(3,1,1);
plot(1:100, Pi_1(11:110), 'r-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,100])
ylabel('Player 16', 'FontSize', 15,'FontWeight','bold');   
title('Player Skill vs. Gibbs Iteration (burn-in & thinning)', 'FontSize', 20,'FontWeight','bold')
   
ax2 = subplot(3,1,2);
plot(1:100, Pi_2(11:110), 'g-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,100])
ylabel('Player 69', 'FontSize', 15,'FontWeight','bold');   

ax3 = subplot(3,1,3);
plot(1:100, Pi_3(11:110), 'b-', 'Linewidth', 1)
set(gca,'fontsize',15);
xlim([0,100])
xlabel('Gibbs Iteration', 'FontSize', 15,'FontWeight','bold');
ylabel('Player 107', 'FontSize', 15,'FontWeight','bold');   
