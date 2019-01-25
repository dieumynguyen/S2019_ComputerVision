
close all; clear;

%%%%% Random seed
rng(2);

mytic = tic;

%%%% Global variables
target = imread('Materials/10x10wifi.jpg');

target_1 = target(:,:,1); % Red channel
target_2 = target(:,:,2); % Green channel
target_3 = target(:,:,3); % Blue channel

[target_x, target_y, target_z] = size(target);

pop_size = numel(target);
channel_size = target_x * target_y;

DNA_bits = [0 255];
max_gen = 50001;
mating_factor = 10;
exponential_factor = 2;
breeding_method = 1;
mutation_rate = 0.01;
kill_factor = 0.9;   % Only breed top x% to save time
tolerance_1 = 50;
tolerance_2 = 50;
mutation_range = [0 30]; % Brighten/darken a pixel by an amount
random_mutation_rate = 0.25;

%%%%%%%%%%%%%%% Save data for plotting
generations = [];
max_fitness_over_t = [];
avg_fitness_over_t = [];
genetic_diverity = [];

%%%%%%%%%%%%%%% Genetic Algorithm 
population = buildPopulation(DNA_bits, pop_size, target_x, target_y, target_z);

img_gens = round(linspace(0, max_gen-1, 9));

gen = 0;
img_num = 1;  % Increment when each of 9 subplots is plotted 
while gen < max_gen
        
    %%%%%%%%%%%%%%%%%% Measure fitness
    % Red channel
    population_fitness_1 = calculateFitness(pop_size, population(:,:,1,:), target_1, ...
                           target_x, target_y, tolerance_1, tolerance_2) + 1e-5;
    [max_fitness_1, max_index_1] = max(population_fitness_1);
    mean_fitness_1 = mean(population_fitness_1);
    
    % Green channel
    population_fitness_2 = calculateFitness(pop_size, population(:,:,2,:), target_2, ...
                           target_x, target_y, tolerance_1, tolerance_2) + 1e-5;
    [max_fitness_2, max_index_2] = max(population_fitness_2);
    mean_fitness_2 = mean(population_fitness_2);

    % Blue channel
    population_fitness_3 = calculateFitness(pop_size, population(:,:,3,:), target_3, ...
                           target_x, target_y, tolerance_1, tolerance_2) + 1e-5;
    [max_fitness_3, max_index_3] = max(population_fitness_3);
    mean_fitness_3 = mean(population_fitness_3);
    
%     %%%%%%%%%%%%%%%%%% Plotting data
%     diversity_1 = max_fitness_1 - mean_fitness_1;
%     max_fitness_over_t = [max_fitness_over_t max_fitness_1];
%     avg_fitness_over_t = [avg_fitness_over_t mean_fitness_1];
%     genetic_diverity = [genetic_diverity diversity_1];
%     generations = [generations gen];
%     
%     %%%%%%%%%%%%%%%%%% Subplot max images
%     if any(img_gens == gen) == 1 
%        subplot(3,3,img_num);
%        image = population(:,:,1,max_index_1);
%        imshow(uint8(image));
%        title(['Generation: ', num2str(gen)]);
%        img_num = img_num + 1;
%     end

    %%%%%%%%%%%%%%%%%% Print report every x generations
    report = ['Gen: ', ...
              num2str(gen, '%.3d'), ...
              ', MaxRed: ', ...
              num2str(max_fitness_1, '%.2f'), ...
              ', AvgRed: ', ...
              num2str(mean_fitness_1, '%.2f'), ...
              ', MaxGreen: ', ...
              num2str(max_fitness_2, '%.2f'), ...
              ', AvgGreen: ', ...
              num2str(mean_fitness_2, '%.2f'), ...
              ', MaxBlue: ', ...
              num2str(max_fitness_3, '%.2f'), ...
              ', AvgBlue: ', ...
              num2str(mean_fitness_3, '%.2f')];
    if mod(gen, 10) == 0
        disp(report)
    end
    
%     %%%%%%%%%%%%%%%%%% Stop loop if target found
%     if max_fitness == 1.00 
%         disp('Exact match found. Ending early...')
%         disp(report)
%         break
%     end
    
    %%%%%%%%%%%%%%%%%% Build mating pool
    % Red channel
    [mating_pool_1, raffles_1] = buildMatingPool(population_fitness_1, ...
                             mating_factor, exponential_factor);
    [mating_pool_2, raffles_2] = buildMatingPool(population_fitness_2, ...
                             mating_factor, exponential_factor);
    [mating_pool_3, raffles_3] = buildMatingPool(population_fitness_3, ...
                             mating_factor, exponential_factor);
    
    %%%%%%%%%%%%%%%%%% Breed new population
    new_population = zeros(size(population));
    for member = 1:pop_size
        %%%%%%%%% Select 2 random parents
        limit_raffle = round(length(mating_pool_1) * (1-kill_factor));
        
        % Red channel
        probs_1 = raffles_1 / sum(raffles_1);
        pick_1_1 = randsample(mating_pool_1(1:limit_raffle), 1, true, probs_1(1:limit_raffle));
        pick_2_1 = randsample(mating_pool_1(1:limit_raffle), 1, true, probs_1(1:limit_raffle));
        parent_1_1 = population(:,:,1,mating_pool_1(pick_1_1));
        parent_2_1 = population(:,:,1,mating_pool_1(pick_2_1));       
        
        % Green channel
        probs_2 = raffles_2 / sum(raffles_2);
        pick_1_2 = randsample(mating_pool_2(1:limit_raffle), 1, true, probs_2(1:limit_raffle));
        pick_2_2 = randsample(mating_pool_2(1:limit_raffle), 1, true, probs_2(1:limit_raffle));
        parent_1_2 = population(:,:,2,mating_pool_2(pick_1_2));
        parent_2_2 = population(:,:,2,mating_pool_2(pick_2_2));  
       
        % Blue channel
        probs_3 = raffles_3 / sum(raffles_3);
        pick_1_3 = randsample(mating_pool_3(1:limit_raffle), 1, true, probs_3(1:limit_raffle));
        pick_2_3 = randsample(mating_pool_3(1:limit_raffle), 1, true, probs_3(1:limit_raffle));
        parent_1_3 = population(:,:,3,mating_pool_3(pick_1_3));
        parent_2_3 = population(:,:,3,mating_pool_3(pick_2_3));  
        
        %%%%%%%%% Breed a progeny 
        progeny_red = breed(parent_1_1, parent_2_1, channel_size, breeding_method);
        progeny_green = breed(parent_1_2, parent_2_2, channel_size, breeding_method);
        progeny_blue = breed(parent_1_3, parent_2_3, channel_size, breeding_method);

        %%%%%%%% Mutate progeny
        progeny_mutation_red = causeMutation(DNA_bits, progeny_red, mutation_rate, ...
                           target_x, target_y, target_z, ...
                           mutation_range, random_mutation_rate);                  
        progeny_mutation_green = causeMutation(DNA_bits, progeny_green, mutation_rate, ...
                   target_x, target_y, target_z, ...
                   mutation_range, random_mutation_rate);                       
        progeny_mutation_blue = causeMutation(DNA_bits, progeny_blue, mutation_rate, ...
                   target_x, target_y, target_z, ...
                   mutation_range, random_mutation_rate);

        new_population(:,:,1,member) = progeny_mutation_red;
        new_population(:,:,2,member) = progeny_mutation_green;
        new_population(:,:,3,member) = progeny_mutation_blue;
    end

    population = new_population;   
                         
    gen = gen + 1;
    
end

end_msg = ['Total generations: ', num2str(gen, '%.3d')];
disp(end_msg)

mytoc = toc(mytic);
% 
% sgtitle(['N = ', num2str(max_gen-1), ', size = ', num2str(target_x), ...
%          'x', num2str(target_y), ', elapsedTime = ', num2str(mytoc), ' sec']);

% %%%%%%%%%%%%%%% Plotting & Writing txt
% % Plot maximum and average fitness over generations
% figure
% plot(max_fitness_over_t)
% hold on
% plot(avg_fitness_over_t)
% plot(genetic_diverity)
% title('Maximum fitness, average fitness, genetic diversity over time')
% xlabel('Time')
% ylabel('Fitness/Diversity')
% xlim([0 max_gen+5]);
% ylim([0 1]);
% legend('Maximum fitness', 'Average fitness', 'Genetic diversity', 'best')
% hold off
%      
% % Write data to text file
% fun_1 = @(x) sprintf('%0.4f', x);
% fun_2 = @(x) sprintf('%0.2e', x);
% fun_3 = @(x) sprintf('%0.3d', x);
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
% gens = num2cell(generations);
% gens_cell = cellfun(fun_3, gens, 'UniformOutput',0);
% gen_T = cell2table(gens_cell', 'VariableNames', {'Gen'});
% 
% T = [gen_T max_T avg_T div_T];
% 
% writetable(T, 'color/wifi_red_channel_report.txt', 'Delimiter', '\t')    

%%%%%%%%%%%%%%% Functions

% 1. Build initial population and store in 4-D matrix
% Access each matrix/img by: population(:,:,:,i)
function population = buildPopulation(target, DNA_bits, pop_size, target_x, target_y, target_z)
%     % Limit colors using frequency/density
%     target_tab = tabulate(unique(target));
%     target_val = target_tab(:,1,:);
%     target_prob = target_tab(:,3,:);

    population = randi(DNA_bits, [target_x, target_y, target_z, pop_size]);
end
    
% 2. Calculate fitness 
% Compare each pop member to target pixel by pixel
function population_fitness = calculateFitness(pop_size, population, target,...
                     target_x, target_y, tolerance_1, tolerance_2)
    % Old way to compare pixel by pixel
    population_fitness = sum(population == target, [1,2,3]) / pop_size; 
    
    % New methods to measure fitness, will have 4 fitness values:
    
    % 1. Percent of values that fit within a tolerance range from the target image
    abs_term_1 = abs(population - double(target));
    fitness_1 = sum(abs_term_1 < tolerance_1, [1,2,3]) / pop_size;
    
    % 2. Check avg values around pixels
    % Compare meanFilter version of target and population with same
    % tolerance
    filtered_target = meanFilter(double(target), target_x, target_y);
    filtered_pop = meanFilter(population, target_x, target_y);
    abs_term_2 = abs(filtered_pop - filtered_target);
    fitness_2 = sum(abs_term_2 < tolerance_1, [1,2,3]) / pop_size;
    
    % 3. Check rate of change bw pixels (how quickly the color changes)
    % Rate of change up/down compared to target image
    diff_target_updown = diff(double(target),1,1);
    diff_pop_updown = diff(population,1,1);
    abs_term_3 = abs(diff_pop_updown - diff_target_updown);
    fitness_3 = sum(abs_term_3 < tolerance_2, [1,2,3]) / pop_size;
    
    % Rate of change left/right compared to target image
    diff_target_lr = diff(double(target),1,1);
    diff_pop_lr = diff(population,1,1);
    abs_term_4 = abs(diff_pop_lr - diff_target_lr);
    fitness_4 = sum(abs_term_4 < tolerance_2, [1,2,3]) / pop_size;

    fit_total = sqrt(fitness_1.^2 + fitness_2.^2 + fitness_3.^2 + fitness_4.^2) / 2;
end

% 2.1. Mean filter to check average values around pixels
function filtered_img = meanFilter(img_to_filter, target_x, target_y)
    h = ones(target_x, target_y);
    filtered_img = imfilter(log(img_to_filter), h, 'replicate');
    filtered_img = exp(filtered_img);
    filtered_img = filtered_img .^ (1/numel(h));
end

% 3. Build mating pool
function [mating_pool, raffles] = buildMatingPool(population_fitness, ...
                                  mating_factor, exponential_factor)
    % Apply exponential factor
    population_fitness_exp = population_fitness.^exponential_factor;    
    % Normalize fitness values
    norm_population_fitness = population_fitness_exp / max(population_fitness_exp);
    % Multiply mating factor
    mating_population_fitness = norm_population_fitness * mating_factor;
    % Put members in pool
    [sorted_raffles, sorted_mating_pool] = sort(mating_population_fitness, 'descend');
    raffles = squeeze(sorted_raffles);
    mating_pool = squeeze(sorted_mating_pool);
end

% 4. Breed
function progeny = breed(parent_1, parent_2, channel_size, breeding_method)
    progeny = parent_1;
    if breeding_method == 0 
        midpoint_dna = randi(numel(parent_1));
        progeny(midpoint_dna:end) = parent_2(midpoint_dna:end);   
    elseif breeding_method == 1
        selections_1_idxs = randi([0 1], 1, channel_size);
        selections_2_idxs = 1 - selections_1_idxs;
        selections_2 = find(selections_2_idxs);
        progeny(selections_2) = parent_2(selections_2);
    end   
end

% 5. Mutate
function progeny_mutation = causeMutation(DNA_bits, progeny, mutation_rate, ...
                            target_x, target_y, target_z, ...
                            mutation_range, random_mutation_rate)
    progeny_mutation = progeny;
                        
    random_sign = randsample([-1 0 1], 1, true, [0.5 0 0.5]);
    offset_mutations_range = randi(mutation_range, [target_x, target_y]);
    mutations_offset = progeny_mutation + random_sign*offset_mutations_range;
    mutations_random = randi(DNA_bits, [target_x, target_y]);
    
    % Choose if pixel will be mutated with mutation rate & mutate with
    % offset
    mask_1 = rand([target_x, target_y]);
    mutations_mask = mask_1 <= mutation_rate; 
    progeny_mutation(mutations_mask) = mutations_offset(mutations_mask);
    
    % Choose if chosen pixel will be RANDOMLY mutated with random mutation rate
    % Offset mutations will be replaced with random mutation
    mask_2 = rand([target_x, target_y]);
    combine_mask = mutations_mask .* mask_2;
    random_mutations_mask = combine_mask <= random_mutation_rate & combine_mask > 0;
    progeny_mutation(random_mutations_mask) = mutations_random(random_mutations_mask);
    
end
