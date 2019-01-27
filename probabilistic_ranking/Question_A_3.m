% Answer for Question A
% Plot auto-correlation

[cov_ww,lags] = xcov(samples(1,:),100,'coeff');

plot(lags,cov_ww,'k-','Linewidth',1);
set(gca,'fontsize',15);
xlim([-100,100])
xlabel('Lag', 'FontSize', 15, 'FontWeight','bold');
ylabel('Auto-correlation', 'FontSize', 15, 'FontWeight','bold');   

title('Auto-correlation vs. Lag (Player 1)', 'FontSize', 20,'FontWeight','bold')
