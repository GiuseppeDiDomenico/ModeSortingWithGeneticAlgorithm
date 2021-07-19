function [score] = evaluation(NObjective,population, objectives)
[m, ~] = size(population);
score = zeros(m, NObjective);

for i=1:m
    score(i,:) = feval(objectives, population(i,:));
end
end