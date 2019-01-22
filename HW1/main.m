
%%%%% Random seed
rng(42);

tic
%%%%% Global variables
pop_size = 200;
target_str = 'To be or not to be';
target_len = strlength(target_str);
DNA_bits = [32, 65:90, 97:122];
max_gen = 500;
mating_factor = 10;
breeding_method = 0;
mutation_rate = 0.01;

%%%%% Save data for plotting
generations = [];
max_fitness_over_t = [];
avg_fitness_over_t = [];
genetic_diverity = [];
fittest_strings = {};

%%%%% Genetic algorithm
% Build initial population 
population = buildPopulation(target_len, pop_size, DNA_bits);

gen = 0;
while gen < max_gen
    % Measure fitness
    population_fitness = calculateFitness(target_str, target_len, population);
    
    % Save data to arrays
    [max_fitness, max_index] = max(population_fitness);
    max_string = population{max_index};
    mean_fitness = mean(population_fitness);
    diversity = max_fitness - mean_fitness;
    max_fitness_over_t = [max_fitness_over_t max_fitness];
    avg_fitness_over_t = [avg_fitness_over_t mean_fitness];
    genetic_diverity = [genetic_diverity diversity];
    fittest_strings{gen+1} = max_string;
    generations = [generations gen];
    
    % Print report every 10 generations
    report = ['Generation: ', num2str(gen, '%.3d'), ', Max fitness: ', ...
              num2str(max_fitness, '%.2f'), ', Fittest string: ', max_string];
    if mod(gen, 10) == 0
        disp(report)
    end
    
    % Stop loop if max fitness = 1.0
    if max_fitness == 1
        disp('Exact match found. Ending early...')
        disp(report)
        break
    end

    % Build mating pool (sorted indices of population members)
    [mating_pool, raffles] = buildMatingPool(population_fitness, mating_factor);

    % New population 
    new_population = {};
    for member = 1:pop_size
        % Select 2 random parents and breed 
        limit_raffle = length(mating_pool) * 0.6;
        pick_1 = randsample(mating_pool(1:limit_raffle), 1, true, raffles(1:limit_raffle));
        pick_2 = randsample(mating_pool(1:limit_raffle), 1, true, raffles(1:limit_raffle));
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

end_msg = ['Total generations: ', num2str(gen, '%.3d')];
disp(end_msg)
toc

% Plot maximum and average fitness over generations
plot(max_fitness_over_t)
hold on
plot(avg_fitness_over_t)
plot(genetic_diverity)
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

writetable(T, 'TEST.txt', 'Delimiter', '\t')
     