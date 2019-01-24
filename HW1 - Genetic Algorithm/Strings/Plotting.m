
cutting = 3;
plot(mutation_rate(1:end-cutting), end_maxes(1:end-cutting))
hold on
plot(mutation_rate(1:end-cutting), end_avgs(1:end-cutting))
xlabel('Mutation rates') 
ylabel('Fitness') 
ylim([0 1])
xticks(mutation_rate)
title('Fitness vs Mutation rates')
legend('Maximum fitness', 'Average fitness')