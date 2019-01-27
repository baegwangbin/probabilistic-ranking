% Answer for Question B

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

psi = inline('normpdf(x)./normcdf(x)');
lambda = inline('(normpdf(x)./normcdf(x)).*( (normpdf(x)./normcdf(x)) + x)');

pv = 0.5;            % prior skill variance (prior mean is always 0)

% initialize matrices of skill marginals - means and precisions
Ms = nan(M,1); 
Ps = nan(M,1);

% initialize matrices of game to skill messages - means and precisions
%Mgs = zeros(N,2); 
%Pgs = zeros(N,2);

% allocate matrices of skill to game messages - means and precisions
Msg = nan(N,2); 
Psg = nan(N,2);

iteration = 30;
Mean_mat = zeros(M, iteration);
Cov_mat = zeros(M, iteration);

Mean_16 = zeros(4, iteration);
Cov_16 = zeros(4, iteration);

pv_list = [0.1, 0.5, 1.0, 2.0];
for pv_iter = 1:4
    %pv = pv_list(pv_iter)
    Mgs = zeros(N,2)+pv_list(pv_iter); 
    Pgs = zeros(N,2)+pv_list(pv_iter);

    for iter=1:iteration
      % (1) compute marginal skills 
      for p=1:M
        % precision first because it is needed for the mean update
        Ps(p) = 1/pv + sum(Pgs(G==p)); 
        Ms(p) = sum(Pgs(G==p).*Mgs(G==p))./Ps(p);
      end

      % (2) compute skill to game messages
      % precision first because it is needed for the mean update
      Psg = Ps(G) - Pgs;
      Msg = (Ps(G).*Ms(G) - Pgs.*Mgs)./Psg;

      % (3) compute game to performance messages
      vgt = 1 + sum(1./Psg, 2);
      mgt = Msg(:,1) - Msg(:,2); % player 1 always wins the way we store data

      % (4) approximate the marginal on performance differences
      Mt = mgt + sqrt(vgt).*psi(mgt./sqrt(vgt));
      Pt = 1./( vgt.*( 1-lambda(mgt./sqrt(vgt)) ) );

      % (5) compute performance to game messages
      ptg = Pt - 1./vgt;
      mtg = (Mt.*Pt - mgt./vgt)./ptg;   

      % (6) compute game to skills messages
      Pgs = 1./(1 + repmat(1./ptg,1,2) + 1./Psg(:,[2 1]));
      Mgs = [mtg, -mtg] + Msg(:,[2 1]);

      Mean_mat(:,iter) = Ms;
      Cov_mat(:,iter) = 1./Ps;
      
      Mean_16(pv_iter, iter) = Mean_mat(16,iter);
      Cov_16(pv_iter, iter) = Cov_mat(16,iter);

    end
end

figure(1)
subplot(2,1,1)


plot(Mean_16(1,:),'k-', 'Linewidth', 2);
hold on; plot(Mean_16(2,:),'r-', 'Linewidth', 2);
hold on; plot(Mean_16(3,:),'g-', 'Linewidth', 2);
hold on; plot(Mean_16(4,:),'b-', 'Linewidth', 2);

legend('k=0','k=0.5','k=1.0','k=2.0')

set(gca,'fontsize',15);
ylabel('Mean', 'FontSize', 15,'FontWeight','bold');   
title('Convergence for Different Initialisation', 'FontSize', 20,'FontWeight','bold')
grid on;


%visualise iterations necessary for convergence

subplot(2,1,2)

plot(Cov_16(1,:),'k-', 'Linewidth', 2);
hold on; plot(Cov_16(2,:),'r-', 'Linewidth', 2);
hold on; plot(Cov_16(3,:),'g-', 'Linewidth', 2);
hold on; plot(Cov_16(4,:),'b-', 'Linewidth', 2);

legend('k=0','k=0.5','k=1.0','k=2.0')

set(gca,'fontsize',15);
ylabel('Variance', 'FontSize', 15,'FontWeight','bold');   
xlabel('EP Iteration', 'FontSize', 15,'FontWeight','bold');   
grid on;



