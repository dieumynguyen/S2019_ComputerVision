
plot(max_gen, max_fitness_overall)
hold on
plot(max_gen, avg_fitness_overall)
xlabel('Mating factor') 
ylabel('Fitness') 
ylim([0 1.1])
xticks(max_gen)
title('Fitness vs mating factor')
legend('Maximum fitness', 'Average fitness')