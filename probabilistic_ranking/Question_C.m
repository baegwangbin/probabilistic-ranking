% Answer for Question C

clear all
clc
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

iteration = 100;
Mean_mat = zeros(M, iteration);
Cov_mat = zeros(M, iteration);

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
  
end

ID = [16,1,5,11];   %Djokovic, Nadal, Federer, Murray

% probability of a player winning against another player
winning_prob = zeros(4,4);
for g=1:4
    for h=1:4
        if g ~= h %Only return if they are different 
            winning_prob(g,h)= normcdf (  (Mean_mat(ID(g),iteration)-Mean_mat(ID(h),iteration)) / sqrt(1 + Cov_mat(ID(g),iteration) + Cov_mat(ID(h),iteration))  ) ;
        end
    end
end

% probability of a player having higher skill than another player
skill_prob = zeros(4,4);
for g = 1:4
    for h = 1:4
        skill_prob(g,h) = 1 - normcdf(0, Mean_mat(ID(g),iteration)-Mean_mat(ID(h),iteration), sqrt(Cov_mat(ID(g),iteration)+Cov_mat(ID(h),iteration)));
    end
end

