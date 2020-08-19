function naming_game_var(w1,a1,num_species,num_flies,mutate_prob,max_iter,seq_length)
% naming game with variable sequence length (min length = 2 bits)
% w1: predation risk weight
% a1: similarity weight
% num_species: number of species (default = 6)
% num_flies: number of flies to simulate per species (default = 10)
% mutate_prob: mutation probability (default = 0.1)
% max_iter: max number of epochs to simulate (default = 15000)

if nargin < 3
    num_species = 6;
end
if nargin < 4
    num_flies = 10;
end
if nargin < 5
    mutate_prob = 0.05;
end
if nargin < 6
    max_iter = 10000;
end
if nargin < 7
    seq_length = 5;
end

min_length = 2;

% initial conditions: all sequences are 1's
sequence_matrix = cell(num_flies,num_species);
for i = 1:num_flies
    for j = 1:num_species
        sequence_matrix{i,j} = ones(1,seq_length);
    end
end
scores = cell(num_flies*num_species,1);

sequences = cell(max_iter);%num_species,max_iter);

fitness_matrix = zeros(num_flies,max_iter,num_species);
fitness_orig = zeros(num_flies,1,num_species);


%%
for t = 1:max_iter
    
    % calculate similarity
    for species = 1:num_species
        for fly = 1:num_flies
            fly1 = sequence_matrix{fly,species};
            for sp2 = (species+1):num_species
                for f2 = 1:num_flies
                    fly2 = sequence_matrix{f2,sp2};
                    len = lcm(length(fly1),length(fly2)); % multiply sequences 
                    fly1 = repmat(fly1,1,len/length(fly1));
                    fly2 = repmat(fly2,1,len/length(fly2));
                    distance = calc_similarity(fly1,fly2);
                    score1 = scores{(species-1)*num_flies + fly};
                    score2 = scores{(sp2-1)*num_flies + f2};
                    score1 = [score1; distance];
                    score2 = [score2; distance];
                    scores{(species-1)*num_flies + fly} = score1;
                    scores{(sp2-1)*num_flies + f2} = score2;
                end
            end
        end
    end
    
    for i = 1:num_species % calculate fitness
        for j = 1:num_flies
            fitness = calc_fitness(sequence_matrix{j,i},w1/(a1+w1));
            score = mean(scores{(i-1)*num_flies+j});
            fitness_matrix(j,t,i) = (a1/(w1+a1))*score + fitness;
        end
    end
    
    
    
    % calculate fitness
    for species = 1:num_species
        scores_temp = zeros(num_flies,1);
        for j = 1:num_flies
            scores_temp(j) = (a1/(a1+w1))*mean(scores{(species-1)*num_flies + j}) + calc_fitness(sequence_matrix{j,species},w1/(a1+w1));
        end
        % find lowest scoring flies and randomly choose one
        ind = find(scores_temp == min(scores_temp));
        ind = ind(randperm(length(ind)));
        % all flies take lowest scoring sequence
        for g = 1:num_flies
            sequence_matrix{g,species} = sequence_matrix{ind(1),species};
        end
        for g = 1:num_flies % clear scores of all flies except the one that was copied (or should it be all?)
%             if g ~= ind(1)
                scores{(species-1)*num_flies+g} = [];
%             end
        end
%         scores{(species-1)*num_flies+(1:num_flies)} = [];
        % mutate
        rr = rand(num_flies,1);
        ri = find(rr < mutate_prob);
        for k = 1:length(ri)
            seq_tm = sequence_matrix{ri(k),species};
            if length(seq_tm) > min_length
                r = randi(4);
                if r == 1
                    seq_tm = transpose_bits(seq_tm);
                elseif r == 2
                    seq_tm = flip_bit(seq_tm);
                elseif r == 3
                    seq_tm = add_bit(seq_tm);
                else
                    seq_tm = delete_bit(seq_tm);
                end
            else
                r = randi(3);
                if r == 1
                    seq_tm = transpose_bits(seq_tm);
                elseif r == 2
                    seq_tm = flip_bit(seq_tm);
                else
                    seq_tm = add_bit(seq_tm);
                end
            end
            sequence_matrix{ri(k),species} = seq_tm;
                
%             if sum(seq_tm) == 0 || sum(seq_tm) == length(seq_tm)
%                 seq_tm = flip_bit(seq_tm);
%             else
%                 r = randi(2);
%                 if r == 1
%                     seq_tm = transpose_bits(seq_tm);
%                 else
%                     seq_tm = flip_bit(seq_tm);
%                 end
%             end
%             sequence_matrix(ri(k),:,species) = seq_tm;
        end
    end
    
    for i = 1:num_species % record sequences
        sequences{t} = sequence_matrix;
    end
    
    
end
%%
filename = ['sequence_w' num2str(w1) '_a' num2str(a1) '_ns' num2str(num_species) '_nf' num2str(num_flies) '_mprob' num2str(mutate_prob)];
fig = figure();
for i = 1:num_species
    sequence1 = mode(sequence_matrix(:,:,i));
    time_array = 0:seq_length;
    plot (time_array,i.*ones(1,seq_length+1),'-k','LineWidth',0.25); hold on;
    flash_on_indexes = find(sequence1==1);
    for flash_i=1:length(flash_on_indexes)
        rectangle('Position',[flash_on_indexes(flash_i)-1,i-0.25,1,0.5],...
            'FaceColor',[1 .9 0],'EdgeColor','k','LineWidth',0.25); hold on;
    end
    
end
ylim([0,num_species+1]);
title(['fitness weight = ' num2str(w1) ', similarity weight = ' num2str(a1) ', mutation probability = ' num2str(mutate_prob)])
set(gca,'plotboxaspectratio',[1,0.5,0.5]);
savefig([filename '.fig']);
print([filename '.png'],'-dpng');
close(fig);

save([filename '.mat']);


end

function longest = calc_similarity(seq1, seq2)
% calculate length of longest common sequence in two wrapped sequences
s1 = horzcat(seq1,seq1(1:end-1));
s2 = horzcat(seq2,seq2(1:end-1));
str_len = zeros(length(s2)-1,1);
for i = 1:length(s2) % shift sequence by one bit
    s2 = circshift(s2,1);
    sdiff = ~(abs(s1 - s2));
    if sum(sdiff) > 0 % if there are some shared bits
        str_len(i) = max(accumarray(nonzeros((cumsum(~sdiff)+1).*sdiff),1));
    else
        str_len(i) = 0;
    end
end
longest = max(str_len)/length(s1);


end

function fitness = calc_fitness(sequence,w1)
% fitness is defined as total time light is on (predation risk)
% PLUS coeff * times the flash switches on (energy cost)
pred_risk = sum(sequence)/length(sequence);
% energy = (numel(find(diff(sequence)==-1))+1)/length(sequence);
% fitness = zeros(size(sequence_matrix,1),1);
% for i = 1:size(sequence_matrix,1)
%     seq = sequence_matrix(i,:);
%     fitness(i) = sum(seq);
% %     robustness = calc_robust(seq);
% end
% w1 = 1;
fitness_model = @(pred_risk) w1 * pred_risk ;
fitness = fitness_model(pred_risk);

end


function sequence = transpose_bits(original_sequence)
% switch two adjacent bits in sequence
ind = randi(length(original_sequence)-1,1);
sequence = original_sequence;
sequence(ind+1) = original_sequence(ind);
sequence(ind) = original_sequence(ind+1);
end

function sequence = flip_bit(original_sequence)
% flip a bit in sequence
ind = randi(length(original_sequence),1);
sequence = original_sequence;
temp = sequence(ind);
temp(:) = ~temp; % changes 0 to 1 and 1 to 0
sequence(ind) = temp;
end

function sequence = add_bit(original_sequence)
ind = randi(length(original_sequence),1);
sequence = horzcat(original_sequence(1:ind),randi(2)-1,original_sequence(ind+1:end));
end

function sequence = delete_bit(original_sequence)
ind = randi(length(original_sequence),1);
sequence = original_sequence;
sequence(ind) = [];
end