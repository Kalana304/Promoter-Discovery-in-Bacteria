// -----------------------------------------------------------------------------------------------------
// BM4321 Genomic Signal Processing
//
// Operations
// 1. Reading of a FASTA file
// 2. Extraction of coding and non-coding DNA
// 3. Reading of a GenBank Protein Table
// 4. Basic analysis of a coding and non-coding region
//
// Objectives
// 1. Familiarization with the coding regions of different organisms (Archaea, Bacteria, Eukaryota etc.)
// 2. Preliminary investigation of gene promoter regions
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2020/05/25
// Kalana Abeywardena, ENTC, University of Moratuwa (160005C@uom.lk) 2020/06/17
// Free to use, distribute and modify for educational purposes with attribution
// -----------------------------------------------------------------------------------------------------

////////////////////Begin of FASTA reading//////////////////////////
function result = remove_eols(text_in)
    // Remove EOLs of fasta file
    keys = find(text_in==10);
    if (isempty(keys)) then
        result = text_in;
    else
    text_out = text_in(1:(keys(1)-1));
    k_n = length(keys)-1;
    for k=1:k_n
        text_i = text_in((keys(k)+1):(keys(k+1)-1));
        text_out = [text_out, text_i];
    end
    result = [text_out,text_in((keys(k_n+1)+1):length(text_in))];
    end
endfunction

function comp=get_comp(base)
    if (base==65) then
        comp = 84;
    elseif (base==67) then
        comp = 71;
    elseif (base==71) then
        comp = 67;
    elseif (base==84) then
        comp = 65;
    end
endfunction

function inc=base_inc(base,c_val)
    if (base==65) then
        inc=c_val+[1 0 0 0]';
    elseif (base==67) then
        inc=c_val+[0 1 0 0]';
    elseif (base==71) then
        inc=c_val+[0 0 1 0]';
    elseif (base==84) then
        inc=c_val+[0 0 0 1]';
    end
endfunction

function gen_code=get_fasta_at(file_name,g_pos,g_end,strand)
   //Estimate the necessary overread to compensate for the EOL charactors of FASTA files
   g_len = g_end-g_pos;
   if g_len>0 then
          n_extra = floor(g_len/70);
          n_offset = floor(g_pos/70);
          file_details = fileinfo(file_name);
          file_len = file_details(1);
          fd = mopen(file_name,'rb');
          mseek(0,fd);
          header = mgetl(fd,1);
          g_start = length(header);
          mseek(g_start+g_pos+n_offset,fd);
          raw_code = mget(g_len+n_extra,'c',fd);
          mclose(fd);
          code_i = remove_eols(raw_code);
          if strand==1 then
              gen_code = code_i;
          else
              //get complementary strand
              len = length(code_i);
              code_c = [];
              for k=1:len
                  code_c = [code_c,get_comp(code_i(k))];
                  gen_code = code_c;
              end
          end
   else
       gen_code = [];
   end

endfunction
/////////////////////END of FASTA reading///////////////////////////

function result=compare_oligo(oligo_a,oligo_b)
    if (length(oligo_a)~=length(oligo_b)) then
        result = -1;
    else
        s = sum(abs(oligo_a-oligo_b))
        if (s==0) then
            result = 1;
        else
            result = 0;
        end
    end
endfunction

////////////////Begin of Protein Table reading//////////////////////
function [gene_array_p, gene_array_n, noncoding_array_p, noncoding_array_n]=get_protein_pos_array(filename)
    fd = csvRead(filename,[],[],"string"); // Reading the protein table in csv format
    [n_entries, _] = size(fd) 
        
    ga_p = [];
    ga_n = [];
    
    nca_p = [];
    nca_n = [];
    
    nc_prev_p = 0; 
    nc_prev_n = 0;
    
    for row=2:n_entries
        p_data = fd(row,:);
        pg_start = strtod(p_data(3));
        pg_stop = strtod(p_data(4));
        
        if (~isempty(p_data(5))) then
            if (strcmp(p_data(5),'-')==0) then
                // disp(sprintf('row (%d) strand:%s ',row,p_data(5)));
                ga_n = [ga_n; pg_start, pg_stop];
                nca_n = [nca_n; nc_prev_n, (pg_start-1)];
                nc_prev_n = pg_stop+1;
            else
                // disp(sprintf('row (%d) strand:%s ',row,p_data(5)));
                ga_p = [ga_p; pg_start, pg_stop];
                nca_p = [nca_p; nc_prev_p, (pg_start-1)];
                nc_prev_p = pg_stop+1;
            end
        end
    end
    
    gene_array_p=ga_p;
    noncoding_array_p=nca_p;
    gene_array_n=ga_n;
    noncoding_array_n=nca_n;
    
endfunction
/////////////////End of Protein Table reading///////////////////////

/////////////////////Begin of FASTA saving//////////////////////////
function save_fasta(filename,header_line,gen_code)
    // Save results into a FASTA file
    fd = mopen(filename,'a');
    mputl('',fd);
    mputl(header_line,fd);
    gen_len = length(gen_code);
    for key=1:gen_len
        mput(gen_code(key),'c');
        if pmodulo(key,70)==0 then
            mputl('',fd);
        end
    end
    mclose(fd)
endfunction
/////////////////////End of FASTA saving////////////////////////////

function save_text(filename,seq_name,data)
    // Save results into a FASTA file
    fd = mopen(filename,'a');
    mfprintf(fd,'%s    %.1f\t%.1f\t%.1f\t%.1f\t%.1f.\n',seq_name,data);
    mclose(fd)
endfunction

function bk=get_base_key(base)
    // Base key as A=1, C=2, G=3, T=4 (alphabetic)
    if (base==65) then
        bk=1;
    elseif (base==67) then
        bk=2;
    elseif (base==71) then
        bk=3;
    elseif (base==84) then
        bk=4;
    // For multiple possibilities assign one for consistency
    // Evenly distributed as much as possible A=3,C=3, G=3, T=2
    elseif (base==ascii('R')) then
        bk=1; // Assign A (A or G)    
    elseif (base==ascii('Y')) then
        bk=2; // Assign C (C or T) 
    elseif (base==ascii('S')) then
        bk=3; // Assign G (C or G) 
    elseif (base==ascii('W')) then
        bk=4; // Assign T (A or T) 
    elseif (base==ascii('K')) then
        bk=3; // Assign G (G or T) 
    elseif (base==ascii('M')) then
        bk=1; // Assign A (A or C) 
    elseif (base==ascii('B')) then
        bk=2; // Assign C (C or G or T) 
    elseif (base==ascii('D')) then
        bk=3; // Assign G (A or G or T) 
    elseif (base==ascii('H')) then
        bk=4; // Assign T (A or C or T) 
    elseif (base==ascii('V')) then
        bk=1; // Assign A (A or C or G)
    elseif (base==ascii('N')) then
        bk=2; // Assign C (any base)
    end
endfunction

// update the position frequency matrix
function pos_freq_mat_out = update_pos_freq_mat(seq,pos_freq_mat_in,set_stop)
    [m,n]=size(pos_freq_mat_in);
    s_len = length(seq);
    stop_point = min([n s_len set_stop]); // Where to stop
    for k=1:stop_point
        pos_base_key = get_base_key(seq(k));
        pos_freq_mat_in(pos_base_key,k)=pos_freq_mat_in(pos_base_key,k)+1;
    end
    pos_freq_mat_out = pos_freq_mat_in;
endfunction

// getting PPM created
function ppm=get_ppm(pos_freq_matrix)
    // Get the position probability matrix
    [m,n]=size(pos_freq_matrix);
    col_sum = sum(pos_freq_matrix,1);
    ppm_temp = [];
    for k=1:n
        ppm_temp = [ppm_temp,pos_freq_matrix(:,k)/col_sum(k)];
    end
    ppm = ppm_temp;
endfunction

// Statistical alignment of a given seq with the obtained PPM
function [pos_scores]=stat_align(seq,ppm)
    // Do a statistical alignment
    [m,n]=size(ppm);
    score_vec = [];
    seq_len = length(seq);
    
    // Do statistical alignment
    for k=1:(seq_len-n+1)
        temp_sum = 0;
        for i=1:n
            pos_base_key = get_base_key(seq(k+i-1));
            temp_sum = temp_sum+log(ppm(pos_base_key,i));
        end
        score_vec = [score_vec,temp_sum];
    end
    pos_scores = score_vec;
endfunction

// Entropy calculation
function [w,su]=ppm_info(ppm,p0)
    [m,n]=size(ppm);
    p_mat = [];
    for k=1:m
        p_row = ones(1,n)/p0(k);
        p_mat = [p_mat; p_row];
    end
    w = log(ppm.*p_mat)/log(2);
    su = ppm.*w;
endfunction

// Statistical alignment of seq with adjusted PPM based on entropy
function [pos_scores,col_entropy,col_keys]=stat_align_entropy(seq,ppm,entropy_thresh)
    // Do a statistical alignment using positions with high entropy
    [m,n]=size(ppm);
    score_vec = [];
    seq_len = length(seq);
    
    // Find entropy of PPM for equiprobable baseline and get columns above the threshold
    [w,su]=ppm_info(ppm,[0.25 0.25 0.25 0.25]);
    col_entropy = sum(su,1); // Get the entropy for each column
    col_pos = 1:n; // Column numbers
    col_keys = col_pos(col_entropy>entropy_thresh); // Get keys of columns above the required entropy threshold
    
    col_key_count = length(col_keys);
    col_key_max = max(col_keys); // Ignore all columns of the PPM after the last column with significant entropy
    
    // Find entropy of consensus to normalize
    con_score = 0;
    for i=1:col_key_count
        con_score = con_score+log(max(ppm(:,col_keys(i)))); // Sum maximum entropy of each significant column
    end
    
    // Do statistical alignment
    for k=1:(seq_len-col_key_max+1)
        temp_sum = 0;
        for i=1:col_key_count
            pos_base_key = get_base_key(seq(k+col_keys(i)-1));
            temp_sum = temp_sum+log(ppm(pos_base_key,col_keys(i)));
        end
        score_vec = [score_vec,temp_sum];
    end
    pos_scores = score_vec;
endfunction

function print_latex_table(num_mat)
// Function to print a LaTeX table
    [m,n]=size(num_mat);
    for i=1:m
        s_out = '';
        for j=1:n
            s_out = sprintf('%s&%1.2f',s_out,num_mat(i,j));
        end
        s_out = sprintf('%s\\\\\\hline',s_out);
        disp(s_out);
    end
endfunction

function [conv_seq] = ppm_to_seq(index_seq)
    out_seq = "";
    i_len = length(index_seq);
    for id=1:i_len
        if (index_seq(id)==1) then
            out_seq = out_seq + "A";
        elseif (index_seq(id)==2) then
            out_seq = out_seq + "C";
        elseif (index_seq(id)==3) then
            out_seq = out_seq + "G";
        else
            out_seq = out_seq + "T";
        end
    end
    conv_seq = out_seq;
endfunction
