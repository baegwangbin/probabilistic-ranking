clear all
clc
load tennis_data

randn('seed',27);         % set the pseudo-random number generator seed

M = size(W,1);            % 107, number of players
N = size(G,1);            % 1801, number of games in 2011 season 

pv = 0.5*ones(M,1);       % prior skill variance 
w = zeros(M,1);           % set skills to prior mean

iteration = 100;
interval = 10;

samples = zeros(M,iteration);                    % raw samples
ind_samples = zeros(M,(iteration/interval));     % sub-sampled

id_1 = 
id_2 = 
id_3 = 

P_1 = zeros(iteration,1);
P_2 = zeros(iteration,1);
P_3 = zeros(iteration,1);

Pi_1 = zeros((iteration/interval),1);
Pi_2 = zeros((iteration/interval),1);
Pi_3 = zeros((iteration/interval),1);


for iter = 1:iteration

  % First, sample performance differences given the skills and outcomes
  
  t = nan(N,1);                     % contains a t_g variable for each game
  for g = 1:N                       % loop over games
    s = w(G(g,1))-w(G(g,2));        % difference in skills
    t(g) = randn()+s;               % performace difference sample
    while t(g) < 0                  % rejection sampling: only positive perf diffs accepted
      t(g) = randn()+s;             % if rejected, sample again
    end
  end 
 
  
  % Second, jointly sample skills given the performance differences
  
  m = nan(M,1);                     % container for the mean of the conditional 
                                    % skill distribution given the t_g samples
  for p = 1:M                       % loop over players
    m(p) = t'*((g==G(:,1)) - (g==G(:,2)));      % (***TO DO***) complete this line
  end
  
  iS = zeros(M,M);                  % container for the sum of precision matrices contributed
                                    % by all the games (likelihood terms)
  for i = 1:M
      for j = 1:i
          if i==j
              iS(i,j) = sum(i==G(:,1)) + sum(i==G(:,2));
          else
              iS(i,j) = -sum((i==G(:,1)).*(j==G(:,2))+(i==G(:,2)).*(j==G(:,1)));
              iS(j,i) = iS(i,j);
          end
      end
  end
      
  iSS = diag(1./pv) + iS; % posterior precision matrix
  % prepare to sample from a multivariate Gaussian
  % Note: inv(M)*z = R\(R'\z) where R = chol(M);
  iR = chol(iSS);  % Cholesky decomposition of the posterior precision matrix
  mu = iR\(iR'\m); % equivalent to inv(iSS)*m but more efficient
    
  % sample from N(mu, inv(iSS))
  w = mu + iR\randn(M,1);

  P_1(iter) = w(id_1);
  P_2(iter) = w(id_2);
  P_3(iter) = w(id_3);
  samples(:,iter) = w;
   
  if mod(iter,10)==0
      Pi_1(iter/10) = w(id_1);
      Pi_2(iter/10) = w(id_2);
      Pi_3(iter/10) = w(id_3);
      ind_samples(:,iter/10) = w;
  end
    
end


figure
ax1 = subplot(3,1,1);
plot(1:100, Pi_20(11:110), 'r-')
set(gca,'fontsize',17);
xlim([0,100])
ylabel('Player 20', 'FontSize', 17,'FontWeight','bold');   
%title('Player Skill vs. Gibbs Iteration', 'FontSize', 17,'FontWeight','bold')
   
ax2 = subplot(3,1,2);
plot(1:100, Pi_40(11:110), 'g-')
set(gca,'fontsize',17);
xlim([0,100])
ylabel('Player 40', 'FontSize', 17,'FontWeight','bold');   
%title('Player 40', 'FontSize', 17,'FontWeight','bold')

ax3 = subplot(3,1,3);
plot(1:100, Pi_60(11:110), 'b-')
set(gca,'fontsize',17);
xlim([0,100])
xlabel('Gibbs Iteration', 'FontSize', 17,'FontWeight','bold');
ylabel('Player 60', 'FontSize', 17,'FontWeight','bold');   
%title('Player 60', 'FontSize', 17,'FontWeight','bold')


%covariance = xcov(samples, 100, 'coeff')