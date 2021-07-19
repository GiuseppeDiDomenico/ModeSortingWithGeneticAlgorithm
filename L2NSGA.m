function [xBest, fval] = L2NSGA(NObjective, popsize, nVariables, fun, iterations, mutationProb, nFrequency, nFOS, initialPop)
% This function implement Linkage Tree based NSGA-II (multi-objective genetic algorithm)
% INPUT
% - NObjective   : number of objectives to optimize
% - popsize      : size of the population
% - nVariables   : size of the chromosomes
% - fun          : objectives functions
% - iterations   : number of generation
% - mutationProb : mutation probability
% - nFrequency   : frequency (n. Generations) for updating the Linkage Learning model
% - nFOS         : number of FOSs to consider for the crossover
% - initialPop   : Initial population; if empty, we create a random new
% population

%% generation 1

if isempty(initialPop)
    population = randi(2, popsize, nVariables);
else 
    population = initialPop;
end

score = evaluation(NObjective, population, fun);
f = non_domination_sort_mod([population, score], NObjective, nVariables);
% In f, genotype, phenotybe, rank and crowding distance are appended as columns
[~, n] = size(f);

%% let's learn the model (linkage tree)
population = f(f(:, n-1)==1, 1:nVariables);
LTModel = linkage(population','average', 'hamming');
[FOS] = extractFOS(LTModel, nFOS);
FOS_ID = unique(FOS); % custers ID

front = f(f(:, n-1)==1, :);
front = front(:, nVariables+1:n);
% h = plot3(front(:,1), front(:,2), front(:, 3), 'bo','MarkerFaceColor', 'b'); %title("Generation 1"), grid
disp('generation 1')
% pause(0.01)

% add results at generation zero
results.generation = 1;
data = f(f(:, n-1)==1, :);
results.xBest = data(:, 1:nVariables);
results.fval = data(:, nVariables+1 : nVariables+NObjective);
save('results.mat', 'results')

%% main loop
for generation=2:iterations
    if mod(generation, nFrequency) == 0
        %% let's learn the model
        population = f(f(:, n-1)==1, 1:nVariables);
        LTModel = linkage(population','average');
        hold off
        [FOS] = extractFOS(LTModel, nFOS);
        FOS_ID = unique(FOS);
    end
    
    if mod(generation, 200) == 0
        results.generation = generation;
        data = f(f(:, n-1)==1, :);
        results.xBest = data(:, 1:nVariables);
        results.fval = data(:, nVariables+1 : nVariables+NObjective);
        save(['results',num2str(generation),'.mat'], 'results')
    end
    
    disp(['generation ' num2str(generation)])
    
    children = zeros(popsize, nVariables);
    [nrow, ~] = size(f);
    
    index = 1;
    for i=1:nrow
        parent = f(i,:);
        donor = tournament_selection(f, 1, 2);
        
        child1 = parent(1, 1:nVariables);
        %bestVal = evaluation(NObjective, parent(1, 1:nNariables), fun);
        
        for k=1:round(nFOS)
            % mask for LT-based crossover
            j =  randi(length(FOS_ID), 1);
            mask = (FOS == FOS_ID(j)); 
            child1(1,mask) = donor(1, mask);
            
            %val = evaluation(NObjective, child1(1, 1:nNariables), fun);
            %if (dominance(val,bestVal) == 1)
            %    child1(1,mask) = parent(1, mask);
            %else
            %    bestVal = val;
            %end
        end
        
        % mutation
        mutationMask = rand(1, nVariables) <= mutationProb; % mutation prob
        child1(:, mutationMask) = ~child1(:, mutationMask);
        
        % add children to the new populaton
        children(index, :) = child1;
        index = index + 1;
    end
    
    f_population = f(:, 1:(n-2)); % ingnore last columns with old rank and crowding distance
    
    % calculate the objectives score for the children
    f_children = [children, evaluation(NObjective, children, fun)];
    
    %% combine old and new population
    union = [f_population; f_children];
    
    %% compute rank dominace + crowding distance
    union = non_domination_sort_mod(union, NObjective, nVariables);
    
    %% select best individual that survive to the next generation
    f = replace_chromosome(union, nVariables, NObjective , popsize);
    
    front = f(f(:, n-1)==1, :); % get non-dominated solution
    front = front(:, (nVariables+1):(nVariables+3)); % get objective scores
    
%    set(h(1),'XData',front(:,1),'YData',front(:,2), 'ZData',front(:,3));
%    title("Generation " + generation)
%    drawnow;
end
front = f(f(:, n-1)==1, :);
xBest = front(:, 1:nVariables);
fval = front(:, (nVariables+1):(nVariables+3));
end

