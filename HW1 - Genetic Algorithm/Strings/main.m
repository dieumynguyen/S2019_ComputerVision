
%%%%% Random seed
rng(42);

%%%% Global variables
pop_size = 200;
% pop_size = 50:50:300;
target_str = 'To be or not to be';
% target_str = 'Meow';
% target_str = 'The universe we observe has precisely the properties we should expect if there is at bottom no design no purpose no evil and no good nothing but blind pitiless indifference DNA neither knows nor cares DNA just is And we dance to its music';
target_len = strlength(target_str);
DNA_bits = [32, 65:90, 97:122];
% % DNA_bits = [32:126];
max_gen = 10000;
mating_factor = 10;
exponential_factor = 3;
breeding_method = 0;
mutation_rate = 0.05;
% mutation_rate = [0, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5];
kill_factor = 0.9;

% %%%%% Genetic algorithm
% end_maxes = [];
% end_avgs = [];
% for i = 1:length(mating_factor)
%     [max_fitness_value, avg_fitness_value] = evolve(pop_size, target_str, ...
%                                              target_len, DNA_bits, max_gen, ...
%                                              mating_factor(i), exponential_factor, ...
%                                              breeding_method, ...
%                                              mutation_rate, kill_factor);
%     end_maxes(i) = max_fitness_value;
%     end_avgs(i) = avg_fitness_value;
% end

% [max_fitness_value, avg_fitness_value] = evolve(pop_size, target_str, ...
%                                              target_len, DNA_bits, max_gen, ...
%                                              mating_factor, exponential_factor, ...
%                                              breeding_method, ...
%                                              mutation_rate, kill_factor);


% function [max_fitness_value, avg_fitness_value] = evolve(pop_size, target_str, ...
%                                                 target_len, DNA_bits, max_gen, ...
%                                                 mating_factor, exponential_factor, ...
%                                                 breeding_method, ...
%                                                 mutation_rate, kill_factor)

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
    if mod(gen, 1000) == 0
        disp(report)
    end

    % Stop loop if max fitness = 1.0 or || strcmp(max_string, target_str) == 1
    if mean_fitness == 1.00
        disp('Exact match found. Ending early...')
        disp(report)
        break
    end

    % Build mating pool (sorted indices of population members)
    [mating_pool, raffles] = buildMatingPool(population_fitness, mating_factor, exponential_factor);

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

disp(max_fitness_value)
disp(avg_fitness_value)

end_msg = ['Total generations: ', num2str(gen, '%.3d')];
disp(end_msg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot maximum and average fitness over generations
plot(max_fitness_over_t)
hold on
plot(avg_fitness_over_t)
plot(genetic_diverity)
title('Maximum fitness, average fitness, genetic diversity over time')
% title('ACSII range: [32:126]');
xlabel('Time')
ylabel('Fitness/Diversity')
ylim([0 1]);
legend('Maximum fitness', 'Average fitness', 'Genetic diversity')
hold off

% % Write data to text file
% fun_1 = @(x) sprintf('%0.4f', x);
% fun_2 = @(x) sprintf('%0.2e', x);
% fun_3 = @(x) sprintf('%0.3d', x);
%
% gens = num2cell(generations);
% gens_cell = cellfun(fun_3, gens, 'UniformOutput',0);
% gen_T = cell2table(gens_cell', 'VariableNames', {'Gen'});
%
% max_fitnesses = num2cell(max_fitness_over_t);
% max_fitnesses_cell = cellfun(fun_1, max_fitnesses, 'UniformOutput',0);
% max_T = cell2table(max_fitnesses_cell', 'VariableNames', {'MaxFit'});
%
% avg_fitnesses = num2cell(avg_fitness_over_t);
% avg_fitnesses_cell = cellfun(fun_1, avg_fitnesses, 'UniformOutput',0);
% avg_T = cell2table(avg_fitnesses_cell', 'VariableNames', {'AvgFit'});
%
% genetic_diverities = num2cell(genetic_diverity);
% genetic_diverities_cell = cellfun(fun_2, genetic_diverities, 'UniformOutput',0);
% div_T = cell2table(genetic_diverities_cell', 'VariableNames', {'Diversity'});
%
% strings = fittest_strings;
% strings_T = cell2table(strings', 'VariableNames', {'FittestString____'});
%
% T = [gen_T strings_T max_T avg_T div_T];
%
% writetable(T, 'FINAL_ORIG_PARAMS_breed1.txt', 'Delimiter', '\t')


% JUNK...
% [max_fitness_overall, avg_fitness_overall] = arrayfun(@(x) evolve(pop_size, target_str, ...
%        target_len, DNA_bits, x, ...
%        mating_factor, breeding_method, mutation_rate, ...
%        kill_factor), max_gen);
% disp(max_fitness_overall)
% disp(avg_fitness_overall)
