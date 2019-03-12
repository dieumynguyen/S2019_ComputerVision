%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1.3 : Build a mating pool
% - Pick 2 parents from the population based on fitness. Fitter 
% phrases have more chance of being picked.
% - Multiple fitness by a mating factor (e.g. 10) and add that 
% many tickets to the raffle
% - Normalize fitness values of each population so the fittest 
% member will always get 10 tickets
% - Function outputs mating pool that contains indices of pop 
% members added to the pool (not the strings)

function [mating_pool, raffles, mating_population_fitness] = buildMatingPool(population_fitness, mating_factor, exponential_factor)
   
    population_fitness = population_fitness.^exponential_factor;    
    
    % Normalize fitness values
    norm_population_fitness = population_fitness / max(population_fitness);
    % Multiply mating factor
    mating_population_fitness = norm_population_fitness * mating_factor;
    % Put members in pool
    [raffles, mating_pool] = sort(mating_population_fitness, 'descend'); 
end