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
