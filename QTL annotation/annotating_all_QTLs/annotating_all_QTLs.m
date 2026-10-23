function annotating_all_QTLs()
    warnState = warning;
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

    %% Read inputs
    X = readtable('final_model_cleaned.csv', 'Delimiter','\t');
    CI = readtable('pos_uncerts.csv');
    A = readtable('essential_genes_Liu_2015.xlsx');
    C = readtable('all_allele_freqs.txt');
    Y=readtable('BYxRM_SNPs_detailed.xlsx'); %list of all SNP positions

    %% Build essential and non-evolvable gene lists
    A_size = size(A);
    essential_gene_list = A{:,1:2};
    non_evolvable_gene_list = strings(1,2);
    count = 1;
    for k = 1:A_size(1)
        if strcmp(A{k,70}, 'Non-evolvable') == 1
            non_evolvable_gene_list(count,1:2) = A{k,1:2};
            count = count + 1;
        end
    end

    %% Filter to additive QTLs only (keep intercept, skip epistatic)
    num_qtls=383; % should be 383 (1 intercept + 382 QTLs)
    X=X(1:num_qtls,:); %only looking at the linear effects for now

    %% Extract lead SNP indices from column 1
    col1 = X{:,1};
    if iscell(col1)
        lead_snp_indices = cellfun(@str2double, col1);
    else
        lead_snp_indices = col1;
    end

    %% Pre-allocate annotation columns
    gene_name_output = strings(num_qtls, 1);
    gene_id_output = strings(num_qtls, 1);
    mutation_type_output = strings(num_qtls, 1);
    nonsense_mutation_state = strings(num_qtls, 1);
    essentiality = strings(num_qtls, 1);
    non_evolvable_essentiality = strings(num_qtls, 1);
    rm_af = repmat(-1, num_qtls, 1);
    maf = repmat(-1, num_qtls, 1);
    intergenic_gene_name_output = strings(num_qtls, 1);
    chr_no_col = repmat(-1, num_qtls, 1);
    chr_pos_col = repmat(-1, num_qtls, 1);
    chr_pos_plotting_col = repmat(-1, num_qtls, 1);
    CI_lead_SNP_col = repmat(-1, num_qtls, 1);
    CI_bottom_SNP_col = repmat(-1, num_qtls, 1);
    CI_top_SNP_col = repmat(-1, num_qtls, 1);
    CI_genes_col = strings(num_qtls, 1);
    CI_genes_with_intergenic_col = strings(num_qtls, 1);

    %% Handle intercept row (row 1, lead SNP = 41594)
    gene_name_output(1) = 'gene_name';
    gene_id_output(1) = 'gene_id';
    mutation_type_output(1) = 'mutation_type';
    nonsense_mutation_state(1) = 'nonsense_mutation_state';
    essentiality(1) = 'gene_essentiality';
    non_evolvable_essentiality(1) = 'non_evolvable_gene_essentiality';
    rm_af(1) = -1;
    maf(1) = -1;
    intergenic_gene_name_output(1) = 'intergenic_gene_name';
    chr_no_col(1) = -1;
    chr_pos_col(1) = -1;
    chr_pos_plotting_col(1) = -1;
    CI_lead_SNP_col(1) = -1;
    CI_bottom_SNP_col(1) = -1;
    CI_top_SNP_col(1) = -1;
    CI_genes_col(1) = 'CI_genes';
    CI_genes_with_intergenic_col(1) = 'CI_genes_with_intergenic';

    %% Annotate each real additive QTL (rows 2 through num_qtls)
    for j = 2:num_qtls
        j
        snp_index_number = lead_snp_indices(j); % 0-indexed lead SNP

        %% Gene name annotation
        [gene_name_output(j), gene_id_output(j)] = gene_name(snp_index_number);

        %% Mutation type annotation
        if strcmp(gene_name_output(j), 'N/A')
            mutation_type_output(j) = 'intergenic';
            nonsense_mutation_state(j) = '0';
        else
            gene_names_split = split(gene_name_output(j), ';');
            gene_ids_split = split(gene_id_output(j), ';');
            mutation_type_temp = '';
            nonsense_temp = '';
            for k = 1:length(gene_names_split)
                [mt, ns] = mutation_type(snp_index_number, 0, gene_names_split{k}, gene_ids_split{k});
                if isempty(mutation_type_temp)
                    mutation_type_temp = mt;
                    nonsense_temp = num2str(ns);
                else
                    mutation_type_temp = append(mutation_type_temp, ';', mt);
                    nonsense_temp = append(nonsense_temp, ';', num2str(ns));
                end
            end
            mutation_type_output(j) = mutation_type_temp;
            nonsense_mutation_state(j) = nonsense_temp;
        end

        %% Essentiality annotation
        gene_names_split = split(gene_name_output(j), ';');
        essentiality_temp = '';
        if strcmp(gene_name_output(j), 'N/A')
            essentiality(j) = 'N/A';
        else
            for k = 1:length(gene_names_split)
                if ismember(gene_names_split(k), essential_gene_list)
                    ess_k = 'essential';
                else
                    ess_k = 'non-essential';
                end
                if isempty(essentiality_temp)
                    essentiality_temp = ess_k;
                else
                    essentiality_temp = append(essentiality_temp, ';', ess_k);
                end
            end
            essentiality(j) = essentiality_temp;
        end

        %% Non-evolvable essentiality annotation
        gene_names_split = split(gene_name_output(j), ';');
        non_evolvable_temp = '';
        if strcmp(gene_name_output(j), 'N/A')
            non_evolvable_essentiality(j) = 'N/A';
        else
            for k = 1:length(gene_names_split)
                if ismember(gene_names_split(k), non_evolvable_gene_list)
                    nev_k = 'non-evolvable essential';
                else
                    nev_k = 'NOT non-evolvable essential';
                end
                if isempty(non_evolvable_temp)
                    non_evolvable_temp = nev_k;
                else
                    non_evolvable_temp = append(non_evolvable_temp, ';', nev_k);
                end
            end
            non_evolvable_essentiality(j) = non_evolvable_temp;
        end

        %% Allele frequency lookup
        snp_id=strcat('snp',num2str(snp_index_number+1));%+1 to go from 0 indexing to the 1 indexing used in the SNP number in BYxRM SNPs and all_allele_freqs
        l=find(strcmp(snp_id,Y{:,2})); 
        corr_curr_QTL_idx=Y{l,1}; %index of qtl in the 0-41609 range indexing used in all_allele_freqs.txt
            
        m=find(C{:,1}==corr_curr_QTL_idx);

        if ~isempty(m)
            rm_af(j) = C{m,6};
            maf(j) = min(C{m,6}, 1 - C{m,6});
        else
            rm_af(j) = -1;
            maf(j) = -1;
        end

        %% Intergenic gene name
        if strcmp(gene_name_output(j), 'N/A') == 1
            intergenic_gene_name_output(j) = intergenic_gene_name(snp_index_number);
        else
            intergenic_gene_name_output(j) = 'genic';
        end

        %% Intronic status check and filtering
        temp = intronic_status(snp_index_number);
        if ~strcmp(temp, 'N/A')
            intronic_gene = extractAfter(temp, 'intronic to ');

            current_gene_names = split(gene_name_output(j), ';');
            is_in_gene_name = ismember(intronic_gene, current_gene_names);

            if is_in_gene_name
                remove_idx = find(strcmp(current_gene_names, intronic_gene));

                current_gene_names(remove_idx) = [];
                if isempty(current_gene_names)
                    gene_name_output(j) = 'N/A';
                else
                    gene_name_output(j) = strjoin(current_gene_names, ';');
                end

                current_gene_ids = split(gene_id_output(j), ';');
                current_gene_ids(remove_idx) = [];
                if isempty(current_gene_ids)
                    gene_id_output(j) = 'N/A';
                else
                    gene_id_output(j) = strjoin(current_gene_ids, ';');
                end

                current_essentiality = split(essentiality(j), ';');
                current_essentiality(remove_idx) = [];
                if isempty(current_essentiality)
                    essentiality(j) = 'N/A';
                else
                    essentiality(j) = strjoin(current_essentiality, ';');
                end

                current_non_evolvable = split(non_evolvable_essentiality(j), ';');
                current_non_evolvable(remove_idx) = [];
                if isempty(current_non_evolvable)
                    non_evolvable_essentiality(j) = 'N/A';
                else
                    non_evolvable_essentiality(j) = strjoin(current_non_evolvable, ';');
                end

                current_mutation_state = split(mutation_type_output(j), ';');
                current_mutation_state(remove_idx) = [];
                if isempty(current_mutation_state)
                    mutation_type_output(j) = 'intergenic';
                else
                    mutation_type_output(j) = strjoin(current_mutation_state, ';');
                end

                current_nonsense_state = split(nonsense_mutation_state(j), ';');
                current_nonsense_state(remove_idx) = [];
                if isempty(current_nonsense_state)
                    nonsense_mutation_state(j) = '0';
                else
                    nonsense_mutation_state(j) = strjoin(current_nonsense_state, ';');
                end
            else
                warning('intronic_status: gene %s from intronic annotation not found in gene_name_output (SNP index %d)', intronic_gene, snp_index_number);
            end

            % append to intergenic_gene_name_output EXCEPT if SNP still belongs to an ORF - keep 'genic', do not append intronic annotation
            if strcmp(intergenic_gene_name_output(j), 'genic') && ~strcmp(gene_name_output(j), 'N/A')
                % SNP still belongs to an ORF - keep 'genic', do not
                % append intronic annotation
            elseif strcmp(intergenic_gene_name_output(j), 'N/A') || strcmp(intergenic_gene_name_output(j), 'genic')
                intergenic_gene_name_output(j) = temp;
            else
                intergenic_gene_name_output(j) = append(intergenic_gene_name_output(j), ';', temp);
            end
        end

        %% Sanity check
        if strcmp(gene_name_output(j), 'N/A') == 1 && strcmp(mutation_type_output(j), 'intergenic') == 0
            fprintf('no gene name is annotated but QTL is also not annotated as intergenic at QTL index j=%d\n', j);
        end

        %% Chromosome position
        [chr_no_col(j), chr_pos_col(j)] = gene_pos_chr_idx(snp_index_number);
        chr_pos_plotting_col(j) = chr_pos_plotting(chr_no_col(j), chr_pos_col(j));

        %% CI columns from pos_uncerts (row j-1 in CI table since CI has no intercept row)
        ci_row = j - 1;
        CI_lead_SNP_col(j) = CI{ci_row, 1};
        CI_bottom_SNP_col(j) = CI{ci_row, 2};
        CI_top_SNP_col(j) = CI{ci_row, 3};

        %% CI gene annotations
        CI_genes_col(j) = CI_genes(CI{ci_row, 2}, CI{ci_row, 3});
        ci_genes_temp = CI_genes_col(j);
        CI_genes_with_intergenic_col(j) = CI_genes_with_intergenic(CI{ci_row, 2}, CI{ci_row, 3}, ci_genes_temp);
    end

    %% Build output table
    X.("gene_name") = gene_name_output;
    X.("gene_id") = gene_id_output;
    X.("mutation_type") = mutation_type_output;
    X.("nonsense_mutation_state") = nonsense_mutation_state;
    X.("gene_essentiality") = essentiality;
    X.("non_evolvable_gene_essentiality") = non_evolvable_essentiality;
    X.("minor_allele_frequency") = maf;
    X.("rm_allele_frequency") = rm_af;
    X.("intergenic_gene_name") = intergenic_gene_name_output;
    X.("chr_no") = chr_no_col;
    X.("chr_pos") = chr_pos_col;
    X.("chr_pos_plotting") = chr_pos_plotting_col;
    X.("CI_lead_SNP") = CI_lead_SNP_col;
    X.("CI_bottom_SNP") = CI_bottom_SNP_col;
    X.("CI_top_SNP") = CI_top_SNP_col;
    X.("CI_genes") = CI_genes_col;
    X.("CI_genes_with_intergenic") = CI_genes_with_intergenic_col;

    %% Write output
    writetable(X, 'output/final_model_cleaned_annotated.txt', 'Delimiter', '\t');
    warning(warnState);
end
