% Answer for Question A
% Run this before running other Question_A_*.m codes

clear all
clc
load tennis_data

randn('seed',27);         % set the pseudo-random number generator seed

M = size(W,1);            % 107, number of players
N = size(G,1);            % 1801, number of games in 2011 season 

pv = 0.5*ones(M,1);       % prior skill variance 
w = zeros(M,1);           % set skills to prior mean

iteration = 1100;
interval = 10;

samples = zeros(M,iteration);                    % raw samples
ind_samples = zeros(M,(iteration/interval));     % sub-sampled

id_1 = 16; %Djokovic
id_2 = 69; %Philipp-Kohlschreiber
id_3 = 107; %Jean-Julien-Rojer

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
    m(p) = t'*((p==G(:,1)) - (p==G(:,2)));      % (***TO DO***) complete this line
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

k = mean(ind_samples,2);    %mean along the rows 
prob = zeros(M,size(ind_samples,2));
Sum = 0;

%Calculate average probaility
for P1=1:M
   for P2=1:M
       for l = 1:size(ind_samples,2)
       prob(P1, l) = normcdf(ind_samples(P1,l) - ind_samples(P2,l));
       end
   end
end

mean_prob= mean(prob,2);

%Plot this ranking system 
[kk, ii] = sort(mean_prob, 'descend');

np = 107;
figure(3)
barh(kk(np:-1:1), 'r')
set(gca,'YTickLabel',W(ii(np:-1:1)),'YTick',1:np,'FontSize',6)
axis([0 1 0.5 np+0.5])
title('Probabilistic Ranking - Average Probability of Winning based on Gibbs Sampling', 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Probability of a player winning against other players', 'FontSize', 12);
ylabel('Player Name', 'FontSize', 12);
