
%%%%% Random seed
rng(42);

%%%% Global variables
pop_size = 200;
target_str = 'To be or not to be';
target_len = strlength(target_str);
DNA_bits = [32, 65:90, 97:122];
max_gen = 300;
mating_factor = 10;
exponential_factor = 3;
breeding_method = 0;
mutation_rate = 0.05;
kill_factor = 0.2;

%%%%% Save data for plotting
generations = [];
max_fitness_over_t = [];
avg_fitness_over_t = [];
genetic_diverity = [];
fittest_strings = {};

% Build initial population
population = buildPopulation(target_len, pop_size, DNA_bits);

gen = 0;
while gen < max_gen
    % Measure fitness
    population_fitness = calculateFitness(target_str, target_len, population);
    population_fitness = population_fitness + 1e-6;


    [max_fitness, max_index] = max(population_fitness);
    max_string = population{max_index};
    mean_fitness = mean(population_fitness);

    diversity = max_fitness - mean_fitness;
    max_fitness_over_t = [max_fitness_over_t max_fitness];
    avg_fitness_over_t = [avg_fitness_over_t mean_fitness];
    genetic_diverity = [genetic_diverity diversity];
    fittest_strings{gen+1} = max_string;
    generations = [generations gen];

    % Print report every x generations
    report = ['Generation: ', num2str(gen, '%.3d'), ...
              ', Max fitness: ', ...
              num2str(max_fitness, '%.2f'), ...
              ', Avg fitness: ', ...
              num2str(mean_fitness, '%.2f'), ...
              ', Fittest string: ', max_string];
    if mod(gen, 10) == 0
        disp(report)
    end

    % Stop loop if target found
    if max_fitness == 1.00
        disp('Exact match found. Ending early...')
        disp(report)
        break
    end

    % Build mating pool (sorted indices of population members)
    [mating_pool, raffles] = buildMatingPool(population_fitness, ...
                             mating_factor, exponential_factor);

    % New population
    new_population = {};
    for member = 1:pop_size
        % Select 2 random parents and breed
        limit_raffle = length(mating_pool) * kill_factor;
        probs = raffles / sum(raffles);
        pick_1 = randsample(mating_pool(1:limit_raffle), 1, true, probs(1:limit_raffle));
        pick_2 = randsample(mating_pool(1:limit_raffle), 1, true, probs(1:limit_raffle));
        % Get corresponding strings
        parent_1 = population{mating_pool(pick_1)};
        parent_2 = population{mating_pool(pick_2)};

        progeny = breed(parent_1, parent_2, target_len, breeding_method);

        % Mutate progeny
        progeny_mutated = causeMutation(progeny, target_len, mutation_rate, DNA_bits);

        % New population with children until all new generation is bred
        new_population{member} = progeny_mutated;
    end

    population = new_population;

    gen = gen + 1;

end

max_fitness_value = max_fitness_over_t(end);
avg_fitness_value = avg_fitness_over_t(end);

end_msg = ['Total generations: ', num2str(gen, '%.3d')];
disp(end_msg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot maximum and average fitness over generations
plot(max_fitness_over_t)
hold on
plot(avg_fitness_over_t)
plot(genetic_diverity)
title('Maximum fitness, average fitness, genetic diversity over time')
xlabel('Time')
ylabel('Fitness/Diversity')
ylim([0 1]);
legend('Maximum fitness', 'Average fitness', 'Genetic diversity')
hold off

% Write data to text file
fun_1 = @(x) sprintf('%0.4f', x);
fun_2 = @(x) sprintf('%0.2e', x);
fun_3 = @(x) sprintf('%0.3d', x);

gens = num2cell(generations);
gens_cell = cellfun(fun_3, gens, 'UniformOutput',0);
gen_T = cell2table(gens_cell', 'VariableNames', {'Gen'});

max_fitnesses = num2cell(max_fitness_over_t);
max_fitnesses_cell = cellfun(fun_1, max_fitnesses, 'UniformOutput',0);
max_T = cell2table(max_fitnesses_cell', 'VariableNames', {'MaxFit'});

avg_fitnesses = num2cell(avg_fitness_over_t);
avg_fitnesses_cell = cellfun(fun_1, avg_fitnesses, 'UniformOutput',0);
avg_T = cell2table(avg_fitnesses_cell', 'VariableNames', {'AvgFit'});

genetic_diverities = num2cell(genetic_diverity);
genetic_diverities_cell = cellfun(fun_2, genetic_diverities, 'UniformOutput',0);
div_T = cell2table(genetic_diverities_cell', 'VariableNames', {'Diversity'});

strings = fittest_strings;
strings_T = cell2table(strings', 'VariableNames', {'FittestString____'});

T = [gen_T strings_T max_T avg_T div_T];

writetable(T, 'FINAL_ORIG_PARAMS_breed1.txt', 'Delimiter', '\t')


%%%%%%%%%%%%%%% Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1.1: Randomly generate initial population of strings
% Target phrase: "To be or not to be"
% - Initialize 200 random strings (same length as target) with 
% letters and spaces (using ASCII and char() typecast)
% - Store the strings in a single cell, which the function outputs
% - Define pop size as a variable

function pop_cell = buildPopulation(target_len, pop_size, DNA_bits)
    pop_cell = {};
    for n = 1:pop_size
        organism_int = randsample(DNA_bits, target_len, true);
        organism_str = char(organism_int);
        pop_cell{n} = organism_str;
    end   
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1.4 Breed a child from 2 parents
% - Combine the DNA (chars) of the 2 parents with 2 methods
% - 1. Choose a random midpoint in 2 phrases, split the 2 parents' DNA
% up accordingly
% - 2. Take random selection of 1 parent's DNA, followed by the remaining
% DNA from the other. 
% - Function take 2 parents and method # (0 or 1), returns a child

function progeny = breed(parent_1, parent_2, target_len, breeding_method)
    if breeding_method == 0
        midpoint_dna = randi([2, target_len - 1]);
        progeny = parent_1;
        progeny(midpoint_dna:end) = parent_2(midpoint_dna:end);
        
    elseif breeding_method == 1
        selections_1_idxs = randi([0 1], 1, 18);
        selections_2_idxs = 1 - selections_1_idxs;
        selections_2 = find(selections_2_idxs);
        progeny = parent_1;
        progeny(selections_2) = parent_2(selections_2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1.5 Cause DNA mutations
% Using mutation rate, mutate the progeny at a dna location

function progeny_mutated = causeMutation(progeny, target_len, mutation_rate, DNA_bits)
    mutate_prob = rand(1);
    progeny_mutated = progeny;

    if mutate_prob <= mutation_rate
        % Random index of progeny to mutate
        mutate_idx = randi(target_len);
        % Random new DNA bit
        new_bit = char(randsample(DNA_bits, 1, true));
        progeny_mutated(mutate_idx) = new_bit;
    end
end

