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

clear;
close();

exec('bm4321_gene_prom_region.sce');
exec('bm4321_sequence_alignment_func.sce');

function [PPM,cons_score]=PPM_reduced_Results(fasta_file,protein_table)

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
    
    n_ppm = 1000; 
    train_seq = [];
    
    for train_key=1:n_ppm
        start_pribnow = pribnow_loc(train_key);
        seq = selected_seq(train_key,:)
        pribnow_row = seq(start_pribnow+1:start_pribnow+pribnow_mat_len);
        //save_fasta('pribnow_train.fasta',sprintf('>sequence %i',train_key),pribnow_row)
        train_seq = [train_seq; pribnow_row];
    end
    
    [train_s,train_b] = size(train_seq);
    
    k_init = 0.01; // Initial probability (to allow log calculations)
    pribnow_pos_freq_matrix = k_init*ones(4,pribnow_mat_len); // Create empty frequency matrix for nucleotides to subsequently generate PPM of pribnow Box
    
    for n_key=1:train_s
        pribnow_post_prom_seq = train_seq(n_key,:);             //disp(ascii(pribnow_post_prom_seq));
        pribnow_pos_freq_matrix = update_pos_freq_mat(pribnow_post_prom_seq,pribnow_pos_freq_matrix,pribnow_mat_len);
    end
    
    prib_ppm = get_ppm(pribnow_pos_freq_matrix); // Generate Pribnow Box PPM
    [prib_b,prib_p]=size(prib_ppm);
    disp("PPM");
    disp(prib_ppm);
    PPM = prib_ppm;
    
    ent_thresh = 0.02;
    
    [pribnow_w,pribnow_su]=ppm_info(prib_ppm,[0.25 0.25 0.25 0.25]);     // Find Pribnow Box entropy
    col_entropy = sum(pribnow_su,1);                                        // Get the entropy for each column
    col_pos = 1:prib_p;                                                          // Column numbers
    col_keys = col_pos(col_entropy>=ent_thresh);                         // Get keys of columns above the required entropy threshold
    red_ppm = prib_ppm(:,col_keys)
    
    [s_ent,idx_ent] = max(red_ppm,'r');
    con_seq_ent =  ppm_to_seq(idx_ent);;
    disp(sprintf("Consensus using reduced PPM : %s",con_seq_ent));
    
    cons_score = sum(log(s_ent));
endfunction



function alignment(fasta_file,protein_table,prib_ppm, psc_ent)
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
    disp(sprintf("Consensus score using reduced PPM : %f",psc_ent));
    
    align_ent = zeros(2,5);
    no_prom_ent = list();
    
    [n_seq,n_pos] = size(selected_seq); disp(sprintf('%s:  : %d',basename(fasta_file),n_seq));
    ent_thresh = 0.02; 
    
    for dec_thresh = -1:-1:-5
        temp_list = [];
        for test_key=1:n_seq
            test_sequence = selected_seq(test_key,:);
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
        align_ent(2,(-1*dec_thresh)) = align_ent(1,(-1*dec_thresh))/n_seq;
        //save_fasta('pribnow_non_aligned_reduced.fasta',sprintf('>threshold: %d',dec_thresh), ascii(string(test_s-align_ent(1,(-1*dec_thresh)))))
    end
    save_text('alignment_results_1.txt',ascii(ascii(fasta_file)(9:21)), align_ent(1,:));
    save_text('alignment_results_1.txt',ascii(ascii(fasta_file)(9:21)), align_ent(2,:)*100);
endfunction

// creating the basis ppm based on the allocated genome sequences
[prib_ppm, psc_ent] = PPM_reduced_Results('Genomes/NZ_AP018808.1.fasta','ProteinTables/NZ_AP018808.1.csv');

// detecting the promoters based on the PPM created earlier
filenames = basename(listfiles('Genomes'));
for i=1:size(filenames)(1)
    fasta_file = 'Genomes/'+ filenames(i) + '.fasta'
    protein_file = 'ProteinTables/'+ filenames(i) + '.csv'
    alignment(fasta_file,protein_file,prib_ppm, psc_ent)
end
