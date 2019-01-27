% Answer for Question A
% Calculate correlation length

cor_len = zeros(107,1);

for i = 1:107
    [cov_ww,lags] = xcov(samples(i,:),100,'coeff');
    cor_len(i) = sum(abs(cov_ww));
end

sum(cor_len)
sum(cor_len)/107
var(cor_len)