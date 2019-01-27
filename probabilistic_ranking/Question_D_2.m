% Answer for Question D
% Probability that a player has higher skill than another player
% Computed from the marginal skills

Mu = zeros(4,1);
Mu(1) = mean(Pi_1);
Mu(2) = mean(Pi_2);
Mu(3) = mean(Pi_3);
Mu(4) = mean(Pi_4);

Var = zeros(4,1);
Var(1) = var(Pi_1);
Var(2) = var(Pi_2);
Var(3) = var(Pi_3);
Var(4) = var(Pi_4);

skill_prob = zeros(4,4);
for g = 1:4
    for h = 1:4
        skill_prob(g,h) = 1 - normcdf(0, Mu(g)-Mu(h), sqrt(Var(g)+Var(h)));
    end
end

skill_prob