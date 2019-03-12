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