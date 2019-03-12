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
