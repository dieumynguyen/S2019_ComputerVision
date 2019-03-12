%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1.2 Calculate fitness of each member in the population
% Fitness score: % of correct chars compared to target phrase

function pop_fitness = calculateFitness(target_str, target_len, population)
    pop_fitness = [];
    % Compare each str in pop_cell against target_str at char level
    for k = 1:length(population)
        organism_str = population{k};
        fitness_score = sum(organism_str == target_str) / target_len;
        pop_fitness(k) = fitness_score;
    end
end
