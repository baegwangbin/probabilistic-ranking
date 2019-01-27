% Answer for Question E
% Ranking - Message Passing

load tennis_data

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

psi = inline('normpdf(x)./normcdf(x)');
lambda = inline('(normpdf(x)./normcdf(x)).*( (normpdf(x)./normcdf(x)) + x)');

pv = 0.5;            % prior skill variance (prior mean is always 0)

% initialize matrices of skill marginals - means and precisions
Ms = nan(M,1); 
Ps = nan(M,1);

% initialize matrices of game to skill messages - means and precisions
Mgs = zeros(N,2); 
Pgs = zeros(N,2);

% allocate matrices of skill to game messages - means and precisions
Msg = nan(N,2); 
Psg = nan(N,2);

no_iterations = 100;
V1 = zeros(M, no_iterations);
V2 = zeros(M, no_iterations);

for iter=1:no_iterations
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
  
  V1(:,iter)=Ms; %Mean matrix per iternation
  V2(:,iter)=1./Ps; %Covariance matrix per iternation
  
end



Mu = V1(:, iter);
V = V2(:, iter);

rank = zeros(M,1);
for i= 1:M
    for j=1:M
        if(not(i==j))
            rank(i) = rank(i) + normcdf( (Mu(i)-Mu(j)) / sqrt(1 + V(i)  + V(j)) ) ;
        end
    end
end

avg_rank = rank./(M-1);

[kk_e, ii_e] = sort(avg_rank, 'descend');
% [kk_v, jj] = sort(vv, 'descend');

np = 107;
figure(3)
barh(kk_e(np:-1:1), 'b')
set(gca,'YTickLabel',W(ii_e(np:-1:1)),'YTick',1:np,'FontSize',5)
axis([0 1 0.5 np+0.5])

title('Message Passing', 'FontSize', 20, 'FontWeight', 'bold')
xlabel('Average Probability of Winning against other Players', 'FontSize', 15, 'FontWeight','bold');
ylabel('Player Names', 'FontSize', 15, 'FontWeight','bold');