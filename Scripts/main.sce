// ------------------------------------------------------------------------------------------- 
// BM4321 Genomic Signal Processing
// By K.G.Abeywardena - 160005C
// Assignment 1 - Promoter Discovery in Bacteria

// 1. Perform a standard local search (for an intact query) to locate a Pribnow box promoter
//    within upstream positions 5 to 30 of each sequence. Using the first 1000 sequences,
//    obtain a position probability matrix (PPM) with 10 positions for the Pribnow box.
// 2. Using a suitable entropy measure, eliminate the redundant positions of (1).
// 3. Perform a statistical alignment for the remaining sequences of (1) using the initial
//    PPM of (1) and the reduced PPM of (2). Compare the two results and hence for the
//    aligned sequences determine the proportion of genes that do not have promoters. For
//    the statistical alignment you may use the thresholds of -1 to -5 (in decrements of 1)
//    w.r.t. the consensus.
// 4. For the genes with no detectable promoter according to (3), perform a local search
//    with a non-intact query to locate any potential mutated promoters. Comment on the
//    result.
// 5. Repeat (1) to (4) for the TTGACA box. This time the sequence should be from 30 to
//    50 positions upstream.
// 6. Using the reduced Pribnow box PPM of (2) perform only the statistical search of (3)
//    for the given thresholds for the remaining 20 genomes of E. coli. Find the proportion
//    of genes that have detectable promoters for each genome.
// --------------------------------------------------------------------------------------------

clc;
clear;
close();


exec('bm4321_gene_prom_region.sce');
exec('bm4321_sequence_alignment_func.sce');

// fasta file and protein table to be read
fasta_file = 'Genomes/Benchmark.fasta';             // change for the genome file required              
protein_table = 'ProteinTables/Benchmark.csv';      // change for the protein file required

[gp,gn,ncp,ncn]= get_protein_pos_array(protein_table);

coding_in_p = size(gp,1); // number of genes in positive strand
coding_in_n = size(gn,1); // number of genes in negative strand

up_thresh_len = 50;     // 50 bases upstream
down_thresh_len = 3;    // 3 bases downstream
safety = 3;             // 3 bases for safety to ensure a 53 bases length sequence is obtained

///////////////////////// Question 01////////////////////////////////

pribnow_query = ascii('TATAAT'); // since performing standard local search
prib_len = 6;

// Get region with Pribnow Box
pribnow_start = 30; // Start position of pribnow Box
pribnow_stop = 5;   // Stop position of pribnow Box

pribnow_len = pribnow_start-pribnow_stop;
pribnow_mat_len = 10; // Length of PPM of pribnow Box

selected_seq = [];
selected_p = [];
selected_n = [];
pribnow_loc = [];

negatives = 0;
non_methionine = 0;
//no_gene_p = 100;

// getting the pribnow sequences from sense strand
for n_key=1:coding_in_p                   // consider the sense strand
    dna_seq = get_fasta_at(fasta_file,gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety-1, 1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len)
    d_len = length(dna_seq); 
    methionine_seq = dna_seq(up_thresh_len+1:up_thresh_len+3); 
    if (methionine_seq == ascii("ATG")) then 
        dna_seq = dna_seq(1:up_thresh_len);
        pribnow_seq = dna_seq((up_thresh_len-pribnow_start+1):(up_thresh_len-pribnow_stop));
        
        if (length(pribnow_seq) == pribnow_len) then 
            [align_x,align_y,prom_pos] = traceback_prom(pribnow_seq,pribnow_query,1,-1,gap_penalty); //disp(sprintf("alignY: %s length %d pos %d",align_y,length(align_y), prom_pos));

            if ((prom_pos>0)&(prom_pos<(pribnow_len-pribnow_mat_len))) then
                dna_row = pribnow_seq(prom_pos+1:prom_pos+pribnow_mat_len);
                if (length(dna_row) == pribnow_mat_len) then
                    //save_fasta('pribnow_sequences.fasta',sprintf('>sequence %i\taligned_position %i',n_key, prom_pos+1),pribnow_seq)
                    selected_seq = [selected_seq; pribnow_seq];
                    selected_p = [selected_p; n_key];
                    pribnow_loc = [pribnow_loc,prom_pos];
                end
            else
                negatives = negatives +1
            end
        end    
    else
        non_methionine = non_methionine + 1;
    end
end

// getting the pribnow sequences from anti-sense strand
for n_key=1:coding_in_n                   // consider the anti-sense strand 
    dna_seq = get_fasta_at(fasta_file,gn(n_key,2)-down_thresh_len+1, gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len)
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq); 
    methionine_seq = dna_seq(up_thresh_len+1:up_thresh_len+3); 
    if (methionine_seq == ascii("ATG")) then 
        dna_seq = dna_seq(1:up_thresh_len);
        pribnow_seq = dna_seq((up_thresh_len-pribnow_start+1):(up_thresh_len-pribnow_stop));
        if (length(pribnow_seq) == pribnow_len) then 
            [align_x,align_y,prom_pos] = traceback_prom(pribnow_seq,pribnow_query,1,-1,gap_penalty);
            if ((prom_pos>0)&(prom_pos<(pribnow_len-pribnow_mat_len))) then
                dna_row = pribnow_seq(prom_pos+1:prom_pos+pribnow_mat_len);
                if (length(dna_row) == pribnow_mat_len) then
                    //save_fasta('pribnow_sequences.fasta',sprintf('>sequence %i\taligned_position %i',n_key+coding_in_p, prom_pos+1),pribnow_seq)
                    selected_seq = [selected_seq; pribnow_seq];
                    selected_n = [selected_n; n_key+coding_in_p];
                    pribnow_loc = [pribnow_loc,prom_pos];
                end
            else
                negatives = negatives +1
            end
        end    
    else
        non_methionine = non_methionine + 1;
    end
end

disp("Question 01");
total_genomes = coding_in_p + coding_in_n;
disp(sprintf("Total genomes available: %i",total_genomes));

total_sel_genomes = length(selected_n)+length(selected_p);
disp(sprintf("Total genomes selected: %i",total_sel_genomes));
disp(sprintf("Total genomes without ATG: %i",non_methionine));

disp(sprintf("Total genomes without Pribnow: %i",negatives));

disp(sprintf("Ratio of selected: %f",(total_sel_genomes/total_genomes)));
disp(sprintf("Ratio of selected from sense strand: %f",(length(selected_p)/coding_in_p)));
disp(sprintf("Ratio of selected from anti-sense strand: %f",(length(selected_n)/coding_in_n)));

total_discarded = negatives + non_methionine; 
disp(sprintf("Ratio of discared: %f",(total_discarded/total_genomes)));
disp(sprintf("Ratio of discared due to no ATG: %f",(non_methionine/total_genomes)));

disp(sprintf("Selected sequence matrix size: %d x %d",(size(selected_seq)(1),size(selected_seq)(2))));

// dividing the selected sequences to train and test
n_ppm = 1000; 
train_seq = [];

for train_key=1:n_ppm
    start_pribnow = pribnow_loc(train_key);
    seq = selected_seq(train_key,:)
    pribnow_row = seq(start_pribnow+1:start_pribnow+pribnow_mat_len);
    //save_fasta('pribnow_train.fasta',sprintf('>sequence %i',train_key),pribnow_row)
    train_seq = [train_seq; pribnow_row];
end

test_seq = selected_seq(n_ppm+1:$,:);

[train_s,train_b] = size(train_seq);
[test_s,test_b] = size(test_seq);
disp(sprintf("Train sequence matrix size: %d x %d",(train_s,train_b)));
disp(sprintf("Test sequence matrix size: %d x %d",(test_s,test_b)));


k_init = 0.01; // Initial probability (to allow log calculations)
pribnow_pos_freq_matrix = k_init*ones(4,pribnow_mat_len); // Create empty frequency matrix for nucleotides to subsequently generate PPM of pribnow Box

for n_key=1:train_s
    pribnow_post_prom_seq = train_seq(n_key,:);             //disp(ascii(pribnow_post_prom_seq));
    pribnow_pos_freq_matrix = update_pos_freq_mat(pribnow_post_prom_seq,pribnow_pos_freq_matrix,pribnow_mat_len);
end

disp("Pribnow Freq. Table");
disp(pribnow_pos_freq_matrix);

prib_ppm = get_ppm(pribnow_pos_freq_matrix); // Generate Pribnow Box PPM
[prib_b,prib_p]=size(prib_ppm);
disp("PPM");
disp(prib_ppm);

/////////////////////////// Question 02//////////////////////////////
ent_thresh = 0.02;

[pribnow_w,pribnow_su]=ppm_info(prib_ppm,[0.25 0.25 0.25 0.25]);     // Find Pribnow Box entropy
col_entropy = sum(pribnow_su,1);                                        // Get the entropy for each column
col_pos = 1:prib_p;                                                          // Column numbers
col_keys = col_pos(col_entropy>=ent_thresh);                         // Get keys of columns above the required entropy threshold
red_ppm = prib_ppm(:,col_keys)

disp("Question 02");

disp("Column entropy values");
disp(col_entropy);

disp("Reduced PPM with entropy")
disp(red_ppm)

////////////////////////// Question 03////////////////////////////////

disp("Question 03");
// Statistical Aligning using Initial PPM

// Gettinng consensus score
[s_ppm,idx_ppm] = max(prib_ppm,'r');
con_seq_ppm = ppm_to_seq(idx_ppm);
disp(sprintf("Consensus using initial PPM : %s",con_seq_ppm));

psc_ppm = sum(log(s_ppm));
disp(sprintf("Consensus score using initial PPM : %f",psc_ppm));

// Aligning
align_ppm = zeros(2,5);
no_prom_ppm = list();

for dec_thresh = -1:-1:-5
    temp_list = [];
    for test_key=1:test_s
        test_sequence = test_seq(test_key,:);
        s_test = stat_align(test_sequence, prib_ppm); //disp(size(s_test));
        r_score = s_test - psc_ppm;
        n_hits = length(r_score(r_score>=dec_thresh));
        if (n_hits > 0) then
            align_ppm(1,(-1*dec_thresh)) = align_ppm(1,(-1*dec_thresh)) + 1;
        else
            temp_list = [temp_list; test_key];
        end
    end
    align_ppm(2,(-1*dec_thresh)) = align_ppm(1,(-1*dec_thresh))/test_s;
    no_prom_ppm(-1*dec_thresh) = temp_list;
    save_fasta('pribnow_non_aligned_initial.fasta',sprintf('>threshold: %d',dec_thresh), ascii(string(test_s-align_ppm(1,(-1*dec_thresh)))))
end

disp("Aligning results with Initial PPM")
disp(align_ppm);

// Statistical Aligning using Reduced PPM

// Gettinng consensus score
[s_ent,idx_ent] = max(red_ppm,'r');
con_seq_ent =  ppm_to_seq(idx_ent);;

disp(sprintf("Consensus using reduced PPM : %s",con_seq_ent));

psc_ent = sum(log(s_ent));
disp(sprintf("Consensus score using reduced PPM : %f",psc_ent));

align_ent = zeros(2,5);
no_prom_ent = list();

for dec_thresh = -1:-1:-5
    temp_list = [];
    for test_key=1:test_s
        test_sequence = test_seq(test_key,:);
        s_test = stat_align_entropy(test_sequence, prib_ppm, ent_thresh)
        r_score = s_test - psc_ent;
        n_hits = length(r_score(r_score>=dec_thresh));
        if (n_hits > 0) then
            align_ent(1,(-1*dec_thresh)) = align_ent(1,(-1*dec_thresh)) + 1;
        else
            temp_list = [temp_list; test_key];
        end
    end
    no_prom_ent(-1*dec_thresh) = temp_list;
    align_ent(2,(-1*dec_thresh)) = align_ent(1,(-1*dec_thresh))/test_s;
    save_fasta('pribnow_non_aligned_reduced.fasta',sprintf('>threshold: %d',dec_thresh), ascii(string(test_s-align_ent(1,(-1*dec_thresh)))))
end
save_text('alignment_results.txt',ascii(ascii(fasta_file)(9:21)), align_ent(1,:));
save_text('alignment_results.txt',ascii(ascii(fasta_file)(9:21)), align_ent(2,:)*100);
//save_text('alignment_results.txt',ascii(ascii(fasta_file)(9:21)), length(test_seq)-align_ent(1,:));
//save_text('alignment_results.txt',ascii(ascii(fasta_file)(9:21)), (1-align_ent(2,:))*100);
disp("Aligning results with Reduced PPM")
disp(align_ent);

////////////////////////// Question 04///////////////////////////////////

disp("Question 04");
n_insertion = 3;
n_mutations = 3;

prib_deletions = zeros(5,5);
prib_matches = list();

for n_mis = 1:5
    prib_matches_i = zeros(5,5);
    mismatch = no_prom_ppm(n_mis);
    l_mis = length(mismatch);
    for m=1:l_mis
        n_test = mismatch(m);
        seq_test = test_seq(n_test,:); 
        l_seq = length(seq_test)
        [align_x,align_y,prom_pos] = traceback_prom(seq_test, pribnow_query,3,-3,gap_penalty); //disp(sprintf("alignY: %s alignX: %s ",align_y,align_x));
        
        if ((prom_pos>0) & (length(align_y)<=l_seq)) then
            prom_len = length(seq_test(prom_pos+1:length(align_y))); 
                        
            if ((l_seq-prom_len+1) > prom_pos) then
                align_x = ascii(align_x);
                test_x = align_x(prom_pos+1:length(align_y)); //Obtain the promoter
                align_y = ascii(align_y);
                test_y = align_y(prom_pos+1:length(align_y)); // obtain the matched promoter
                
                if (length(test_x)>= prib_len)&(length(test_x) == length(test_y))then
                    deletions = length(test_x(test_x == ascii('-')));
                    
                    if (deletions > 0) then
                        if (deletions <= n_insertion) then
                            prib_deletions(n_mis,deletions+1) =  prib_deletions(n_mis,deletions+1) + 1;
                        elseif (deletions > n_insertion) then
                            prib_deletions(n_mis,5) =  prib_deletions(n_mis,5) + 1;
                        end
                    else
                        [n_inst, n_mut] = pribnow_mutate(test_x, test_y); //disp(sprintf("n_inst: %d n_mut: %d ",n_inst, n_mut));
                        
                        if ((n_inst <= n_insertion)&(n_mut <= n_mutations)) then 
                            prib_matches_i(n_inst+1,n_mut+1) = prib_matches_i(n_inst+1,n_mut+1) + 1;
                        elseif ((n_inst > n_insertion)&(n_mut <= n_mutations)) then 
                            prib_matches_i(5,n_mut+1) = prib_matches_i(5,n_mut+1) + 1;
                        elseif ((n_inst <= n_insertion)&(n_mut > n_mutations)) then 
                            prib_matches_i(n_inst+1,5) = prib_matches_i(n_inst+1,5) + 1;
                        else
                            prib_matches_i(5,5) = prib_matches_i(5,5) + 1;
                        end
                     end
                 end  
             end
         end
         prib_matches(n_mis) = prib_matches_i;
     end
end

//save_fasta('pribnow_init_mutated.fasta',sprintf('>Gap_Count: threshold %d',-1), ascii(prib_gaps_init));

prib_del_ent = zeros(5,5);
prib_matches_ent = list();

for n_mis = 1:5
    prib_matches_i = zeros(5,5);
    mismatch = no_prom_ent(n_mis);
    l_mis = length(mismatch);
    for m=1:l_mis
        n_test = mismatch(m);
        seq_test = test_seq(n_test,:); 
        l_seq = length(seq_test)
        [align_x,align_y,prom_pos] = traceback_prom(seq_test, pribnow_query,3,-3,gap_penalty); //disp(sprintf("alignY: %s alignX: %s ",align_y,align_x));
        
        if ((prom_pos>0) & (length(align_y)<=l_seq)) then
            prom_len = length(seq_test(prom_pos+1:length(align_y))); 
                        
            if ((l_seq-prom_len+1) > prom_pos) then
                align_x = ascii(align_x);
                test_x = align_x(prom_pos+1:length(align_y)); //Obtain the promoter
                align_y = ascii(align_y);
                test_y = align_y(prom_pos+1:length(align_y)); // obtain the matched promoter
                
                if (length(test_x)>= prib_len)&(length(test_x) == length(test_y))then
                    deletions = length(test_x(test_x == ascii('-')));
                    
                    if (deletions > 0) then
                        if (deletions <= n_insertion) then
                            prib_del_ent(n_mis,deletions+1) =  prib_del_ent(n_mis,deletions+1) + 1;
                        elseif (deletions > n_insertion) then
                            prib_del_ent(n_mis,5) =  prib_del_ent(n_mis,deletions) + 1;
                        end
                    else
                        [n_inst, n_mut] = pribnow_mutate(test_x, test_y); //disp(sprintf("n_inst: %d n_mut: %d ",n_inst, n_mut));
                        
                        if ((n_inst <= n_insertion)&(n_mut <= n_mutations)) then 
                            prib_matches_i(n_inst+1,n_mut+1) = prib_matches_i(n_inst+1,n_mut+1) + 1;
                        elseif ((n_inst > n_insertion)&(n_mut <= n_mutations)) then 
                            prib_matches_i(5,n_mut+1) = prib_matches_i(5,n_mut+1) + 1;
                        elseif ((n_inst <= n_insertion)&(n_mut > n_mutations)) then 
                            prib_matches_i(n_inst+1,5) = prib_matches_i(n_inst+1,5) + 1;
                        else
                            prib_matches_i(5,5) = prib_matches_i(5,5) + 1;
                        end
                     end
                 end  
             end
         end
         prib_matches_ent(n_mis) = prib_matches_i;
     end
end

prib_hits = zeros(6,5);

for th = 1:5
    cont_table = prib_matches(th);
    total = sum(cont_table);
    hits = sum(cont_table(1:3,1:3))-cont_table(1,1);
    prib_hits(1,th) = total;
    prib_hits(2,th) = hits;
    prib_hits(3,th) = hits/total;
end

for th = 1:5
    cont_table = prib_matches_ent(th);
    total = sum(cont_table);
    hits = sum(cont_table(1:3,1:3))-cont_table(1,1);
    prib_hits(4,th) = total;
    prib_hits(5,th) = hits;
    prib_hits(6,th) = hits/total;
end

disp("Mutated Pribnow detection - Initial PPM")
disp(prib_hits(2:3,:));

disp("Mutated Pribnow detection - Reduced PPM")
disp(prib_hits(5:6,:));

save_text('mutated_results.txt','Total detections: Initial', prib_hits(1,:));
save_text('mutated_results.txt','Mutated detections: Initial', prib_hits(2,:));
save_text('mutated_results.txt','Percentage detecttions: Initial', prib_hits(3,:)*100);
save_text('mutated_results.txt','Total detections: Reduced', prib_hits(4,:));
save_text('mutated_results.txt','Mutated detections: Reduced', prib_hits(5,:));
save_text('mutated_results.txt','Percentage detecttions: Reduced', prib_hits(6,:)*100);


/////////////////////////// Question 05/////////////////////////////

sigma_query = ascii("TTGACA");
sig_len = 6;

sigma_start = 50;
sigma_stop = 30;

sigma_len = sigma_start-sigma_stop;
sigma_mat_len = 10; // Length of PPM of Sigma Factor

selected_seq = [];
selected_p = [];
selected_n = [];
sigma_loc = [];

negatives = 0;
non_methionine = 0;
//no_gene_p = 100;

// getting the pribnow sequences from sense strand
for n_key=1:coding_in_p                   // consider the sense strand
    dna_seq = get_fasta_at(fasta_file,gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety-1, 1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len)
    d_len = length(dna_seq); 
    methionine_seq = dna_seq(up_thresh_len+1:up_thresh_len+3); 
    if (methionine_seq == ascii("ATG")) then 
        dna_seq = dna_seq(1:up_thresh_len);
        sigma_seq = dna_seq((up_thresh_len-sigma_start+1):(up_thresh_len-sigma_stop));
        
        if (length(sigma_seq) == sigma_len) then 
            [align_x,align_y,sigma_pos] = traceback_local(sigma_seq, sigma_query,1,-1,gap_penalty); //disp(sprintf("alignY: %s length %d pos %d",align_y,length(align_y), prom_pos));
            
            if ((sigma_pos>0)&(sigma_pos<(sigma_len-sigma_mat_len))) then
                dna_row = sigma_seq(sigma_pos+1:sigma_pos+sigma_mat_len);
                if (length(dna_row) == sigma_mat_len) then
                    save_fasta('sigma_sequences.fasta',sprintf('>sequence %i\taligned_position %i',n_key, sigma_pos+1),sigma_seq);//disp(sprintf("alignY: %s alignX %s ",align_y, align_x));
                    selected_seq = [selected_seq; sigma_seq];
                    selected_p = [selected_p; n_key];
                    sigma_loc = [sigma_loc,sigma_pos];
                end
            else
                negatives = negatives +1
            end
        end    
    else
        non_methionine = non_methionine + 1;
    end
end

// getting the pribnow sequences from anti-sense strand
for n_key=1:coding_in_n                   // consider the anti-sense strand 
    dna_seq = get_fasta_at(fasta_file,gn(n_key,2)-down_thresh_len+1, gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len)
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq); 
    methionine_seq = dna_seq(up_thresh_len+1:up_thresh_len+3); 
    if (methionine_seq == ascii("ATG")) then 
        dna_seq = dna_seq(1:up_thresh_len);
        sigma_seq = dna_seq((up_thresh_len-sigma_start+1):(up_thresh_len-sigma_stop));
        
        if (length(sigma_seq) == sigma_len) then 
            [align_x,align_y,sigma_pos] = traceback_local(sigma_seq, sigma_query,1,-1,gap_penalty); //disp(sprintf("alignY: %s length %d pos %d",align_y,length(align_y), prom_pos));
            
            if ((sigma_pos>0)&(sigma_pos<(sigma_len-sigma_mat_len))) then
                dna_row = sigma_seq(sigma_pos+1:sigma_pos+sigma_mat_len);
                if (length(dna_row) == sigma_mat_len) then
                    save_fasta('sigma_sequences.fasta',sprintf('>sequence %i\taligned_position %i',n_key+coding_in_p, sigma_pos+1),sigma_seq)
                    selected_seq = [selected_seq; sigma_seq];
                    selected_n = [selected_n; n_key+coding_in_p];
                    sigma_loc = [sigma_loc,sigma_pos];
                end
            else
                negatives = negatives +1
            end
        end    
    else
        non_methionine = non_methionine + 1;
    end
end

disp("Question 05");
total_genomes = coding_in_p + coding_in_n;
disp(sprintf("Total genomes available: %i",total_genomes));

total_sel_genomes = length(selected_n)+length(selected_p);
disp(sprintf("Total genomes selected: %i",total_sel_genomes));

disp(sprintf("Ratio of selected: %f",(total_sel_genomes/total_genomes)));
disp(sprintf("Ratio of selected from sense strand: %f",(length(selected_p)/coding_in_p)));
disp(sprintf("Ratio of selected from anti-sense strand: %f",(length(selected_n)/coding_in_n)));

total_discarded = negatives + non_methionine; 
disp(sprintf("Ratio of discared: %f",(total_discarded/total_genomes)));
disp(sprintf("Ratio of discared due to no ATG: %f",(non_methionine/total_genomes)));

disp(sprintf("Selected sequence matrix size: %d x %d",(size(selected_seq)(1),size(selected_seq)(2))));

// dividing the selected sequences to train and test
n_ppm = 1000; 
train_seq = [];

for train_key=1:n_ppm
    start_sigma = sigma_loc(train_key);
    seq = selected_seq(train_key,:)
    sigma_row = seq(start_sigma+1:start_sigma+sigma_mat_len);
    save_fasta('sigma_train.fasta',sprintf('>sequence %i',train_key),sigma_row)
    train_seq = [train_seq; sigma_row];
end

test_seq = selected_seq(n_ppm+1:$,:);

[train_s,train_b] = size(train_seq);
[test_s,test_b] = size(test_seq);
disp(sprintf("Train sequence matrix size: %d x %d",(train_s,train_b)));
disp(sprintf("Test sequence matrix size: %d x %d",(test_s,test_b)));


k_init = 0.01; // Initial probability (to allow log calculations)
sigma_pos_freq_matrix = k_init*ones(4,sigma_mat_len); // Create empty frequency matrix for nucleotides to subsequently generate PPM of pribnow Box

for n_key=1:train_s
    sigma_post_prom_seq = train_seq(n_key,:);             //disp(ascii(pribnow_post_prom_seq));
    sigma_pos_freq_matrix = update_pos_freq_mat(sigma_post_prom_seq,sigma_pos_freq_matrix,sigma_mat_len);
end

disp("Sigma Freq. Table");
disp(sigma_pos_freq_matrix);

sigma_ppm = get_ppm(sigma_pos_freq_matrix); // Generate Pribnow Box PPM
[sigma_b,sigma_p]=size(sigma_ppm);
disp("PPM");
disp(sigma_ppm);

ent_thresh = 0.01;

[sigma_w,sigma_su]=ppm_info(sigma_ppm,[0.25 0.25 0.25 0.25]);     // Find sigma binding site entropy
col_entropy = sum(sigma_su,1);                                      // Get the entropy for each column
col_pos = 1:sigma_p;                                               // Column numbers
col_keys = col_pos(col_entropy>=ent_thresh);                       // Get keys of columns above the required entropy threshold
sigma_red_ppm = sigma_ppm(:,col_keys);

disp("Column entropy values");
disp(col_entropy);

disp("Reduced PPM with entropy");
disp(sigma_red_ppm);

// Statistical Aligning using Initial PPM

// Gettinng consensus score
[s_ppm,idx_ppm] = max(sigma_ppm,'r');
con_seq_ppm = ppm_to_seq(idx_ppm);
disp(sprintf("Consensus using initial PPM : %s",con_seq_ppm));

psc_ppm = sum(log(s_ppm));
disp(sprintf("Consensus score using initial PPM : %f",psc_ppm));

// Aligning
align_ppm = zeros(2,5);
no_sigma_ppm = list();

for dec_thresh = -1:-1:-5
    temp_list = [];
    for test_key=1:test_s
        test_sequence = test_seq(test_key,:);
        s_test = stat_align(test_sequence, sigma_ppm); //disp(size(s_test));
        r_score = s_test - psc_ppm;
        n_hits = length(r_score(r_score>=dec_thresh));
        if (n_hits > 0) then
            align_ppm(1,(-1*dec_thresh)) = align_ppm(1,(-1*dec_thresh)) + 1;
        else
            temp_list = [temp_list; test_key]; 
        end
    end
    align_ppm(2,(-1*dec_thresh)) = align_ppm(1,(-1*dec_thresh))/test_s;
    no_sigma_ppm(-1*dec_thresh) = temp_list;
    save_fasta('sigma_non_aligned_initial.fasta',sprintf('>threshold: %d',dec_thresh), ascii(string(test_s-align_ppm(1,(-1*dec_thresh)))))
end

disp("Aligning results with Initial PPM")
disp(align_ppm);

// Statistical Aligning using Reduced PPM

// Gettinng consensus score
[s_ent,idx_ent] = max(sigma_red_ppm,'r');
con_seq_ent =  ppm_to_seq(idx_ent);;

disp(sprintf("Consensus using reduced PPM : %s",con_seq_ent));

psc_ent = sum(log(s_ent));
disp(sprintf("Consensus score using reduced PPM : %f",psc_ent));

align_ent = zeros(2,5);
no_sigma_ent = list();

for dec_thresh = -1:-1:-5
    temp_list = [];
    for test_key=1:test_s
        test_sequence = test_seq(test_key,:);
        s_test = stat_align_entropy(test_sequence, sigma_ppm, ent_thresh);//disp(s_test)
        r_score = s_test - psc_ent;//disp(r_score)
        n_hits = length(r_score(r_score>=dec_thresh));
        if (n_hits > 0) then
            align_ent(1,(-1*dec_thresh)) = align_ent(1,(-1*dec_thresh)) + 1;
        else
            temp_list = [temp_list; test_key];
        end
    end
    no_sigma_ent(-1*dec_thresh) = temp_list;
    align_ent(2,(-1*dec_thresh)) = align_ent(1,(-1*dec_thresh))/test_s;
    save_fasta('sigma_non_aligned_reduced.fasta',sprintf('>threshold: %d',dec_thresh), ascii(string(test_s-align_ent(1,(-1*dec_thresh)))))
end

disp("Aligning results with Reduced PPM")
disp(align_ent);

n_insertion = 3;
n_mutations = 3;

sig_deletions = zeros(2,5);
sig_pos_mutations = zeros(1,6);
sig_point_mut = zeros(1,2);
sig_matches = list();

for n_mis = 1:2
    //disp(n_mis);
    sig_matches_i = zeros(5,5);
    mismatch = no_sigma_ppm(n_mis);
    l_mis = length(mismatch);
    mutated_sigma = [];
    for m=1:l_mis
        n_test = mismatch(m);
        seq_test = test_seq(n_test,:); 
        l_seq = length(seq_test)
        [align_x,align_y,sigma_pos] = traceback_local(seq_test, sigma_query,3,-3,gap_penalty); //disp(sprintf("alignY: %s alignX %s",align_y,align_x));
        if (length(align_y)<=l_seq) then
            sigma_len2 = length(seq_test(sigma_pos+1:length(align_y)));
            //promoter_lengths = [promoter_lengths,prom_len2];
            
            if ((l_seq-sigma_len2+1) > sigma_pos) then
                align_x = ascii(align_x);
                test_x = align_x(sigma_pos+1:length(align_y))
                align_y = ascii(align_y);
                test_y = align_y(sigma_pos+1:length(align_y)); //disp(sprintf("testY: %s testX %s",ascii(test_y),ascii(test_x)));// obtain the matched promoter
                
                if (length(test_x)>= sig_len)&(length(test_x) == length(test_y))then
                    deletions = length(test_x(test_x == ascii('-'))); 
                    
                    if (deletions > 0) then
                        if (deletions <= n_insertion) then
                            sig_deletions(n_mis,deletions+1) =  sig_deletions(n_mis,deletions+1) + 1;
                        elseif (deletions > n_insertion) then
                            sig_deletions(n_mis,5) =  sig_deletions(n_mis,5) + 1;
                        end
                    else
                        [pos_arr, point_arr, n_inst, n_mut] = sigma_mutate(test_x, test_y); //disp(sprintf("testY: %s testX %s",ascii(test_y),ascii(test_x)))//disp(sprintf("n_inst: %d n_mut: %d ",n_inst, n_mut));
                        if ((n_inst <= n_insertion)&(n_mut <= n_mutations)) then 
                            sig_matches_i(n_inst+1,n_mut+1) = sig_matches_i(n_inst+1,n_mut+1) + 1;
                        elseif ((n_inst > n_insertion)&(n_mut <= n_mutations)) then 
                            sig_matches_i(5,n_mut+1) = sig_matches_i(5,n_mut+1) + 1;
                        elseif ((n_inst <= n_insertion)&(n_mut > n_mutations)) then 
                            sig_matches_i(n_inst+1,5) = sig_matches_i(n_inst+1,5) + 1;
                        else
                            sig_matches_i(5,5) = sig_matches_i(5,5) + 1;
                        end
                        
                        if (n_mis == 1) then
                            sig_pos_mutations = sig_pos_mutations + pos_arr;
                            sig_point_mut = sig_point_mut + point_arr;
                        end
                     end
                 end  
                
             end
         end
    end
    sig_matches(n_mis) = sig_matches_i;
end

//save_fasta('sigma_init_mutated.fasta',sprintf('>Gap_Count: threshold %d',-1), ascii(sigma_gaps_init));

sig_deletions_ent = zeros(2,5);
sig_pos_mutations_ent = zeros(1,6);
sig_point_mut_ent = zeros(1,2);
sig_matches_ent = list();

for n_mis = 1:2
    sig_matches_i = zeros(5,5);
    mismatch = no_sigma_ent(n_mis);
    l_mis = length(mismatch);
    mutated_sigma = [];
    for m=1:l_mis
        n_test = mismatch(m);
        seq_test = test_seq(n_test,:); 
        l_seq = length(seq_test)
        [align_x,align_y,sigma_pos] = traceback_local(seq_test, sigma_query,3,-3,gap_penalty); //disp(sprintf("alignY: %s alignX %s",align_y,align_x));
        if (length(align_y)<=l_seq) then
            sigma_len2 = length(seq_test(sigma_pos+1:length(align_y)));
            //promoter_lengths = [promoter_lengths,prom_len2];
            
            if ((l_seq-sigma_len2+1) > sigma_pos) then
                align_x = ascii(align_x);
                test_x = align_x(sigma_pos+1:length(align_y))
                align_y = ascii(align_y);
                test_y = align_y(sigma_pos+1:length(align_y)); //disp(sprintf("testY: %s testX %s",ascii(test_y),ascii(test_x)));// obtain the matched promoter
                
                if (length(test_x)>= sig_len)&(length(test_x) == length(test_y))then
                    deletions = length(test_x(test_x == ascii('-'))); //disp(deletions);
                    
                    if (deletions > 0) then
                        if (deletions <= n_insertion) then
                            sig_deletions_ent(n_mis,deletions+1) =  sig_deletions_ent(n_mis,deletions+1) + 1;
                        elseif (deletions > n_insertion) then
                            sig_deletions_ent(n_mis,5) =  sig_deletions_ent(n_mis,5) + 1;
                        end
                    else
                        [pos_arr, point_arr, n_inst, n_mut] = sigma_mutate(test_x, test_y); //disp(sprintf("testY: %s testX %s",ascii(test_y),ascii(test_x)))//disp(sprintf("n_inst: %d n_mut: %d ",n_inst, n_mut));
                        if ((n_inst <= n_insertion)&(n_mut <= n_mutations)) then 
                            sig_matches_i(n_inst+1,n_mut+1) = sig_matches_i(n_inst+1,n_mut+1) + 1;
                        elseif ((n_inst > n_insertion)&(n_mut <= n_mutations)) then 
                            sig_matches_i(5,n_mut+1) = sig_matches_i(5,n_mut+1) + 1;
                        elseif ((n_inst <= n_insertion)&(n_mut > n_mutations)) then 
                            sig_matches_i(n_inst+1,5) = sig_matches_i(n_inst+1,5) + 1;
                        else
                            sig_matches_i(5,5) = sig_matches_i(5,5) + 1;
                        end
                        
                        if (n_mis == 1) then
                            sig_pos_mutations_ent = sig_pos_mutations_ent + pos_arr;
                            sig_point_mut_ent = sig_point_mut_ent + point_arr;
                        end
                     end
                 end  
                
             end
         end
    end
    sig_matches_ent(n_mis) = sig_matches_i;
end

sig_hits = zeros(6,2);

for th = 1:2
    cont_table = sig_matches(th);
    total = sum(cont_table);
    hits = sum(cont_table(1:3,1:3))-cont_table(1,1);
    sig_hits(1,th) = total;
    sig_hits(2,th) = hits;
    sig_hits(3,th) = hits/total;
end

for th = 1:2
    cont_table = sig_matches_ent(th);
    total = sum(cont_table);
    hits = sum(cont_table(1:3,1:3))-cont_table(1,1);
    sig_hits(4,th) = total;
    sig_hits(5,th) = hits;
    sig_hits(6,th) = hits/total;
end

disp("Mutated Pribnow detection - Initial PPM")
disp(sig_hits(2:3,:));

disp("Mutated Pribnow detection - Reduced PPM")
disp(sig_hits(5:6,:));
//save_fasta('sigma_reduced_mutated.fasta',sprintf('>Gap_Count: threshold %d',-1), ascii(sigma_gaps_ent));

