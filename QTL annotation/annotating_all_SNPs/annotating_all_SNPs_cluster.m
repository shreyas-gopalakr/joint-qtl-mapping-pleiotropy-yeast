function annotating_all_SNPs_cluster(i) %input i tells you which SNPs to annotate (let's do i*1000 to (i+1)*1000
    warnState = warning;
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    X=readtable('BYxRM_SNPs_detailed.xlsx'); %list of SNP positions
    A=readtable('essential_genes_Liu_2015.xlsx');
    end_index=height(X);
    step_size=1000; %change step size to 1000 when needed
    starting_snp = (i*step_size); %i should range from 0 to 41
    ending_snp = min(end_index-1, starting_snp + step_size - 1);
    %starting_snp, ending_snp, end_index are all 0-indexed

    A_size=size(A);
    essential_gene_list=A{:,1:2}; %newly changed
    non_evolvable_gene_list=strings(1,2); %newly changed
    count=1;
    for k=1:A_size(1)
        if strcmp(A{k,70},'Non-evolvable')==1
            non_evolvable_gene_list(count,1:2)=A{k,1:2}; %newly changed
            count=count+1;
        end
    end


    C=readtable("all_allele_freqs.txt");

    gene_name_output=strings(ending_snp-starting_snp+1,1);
    gene_id_output=strings(ending_snp-starting_snp+1,1); %newly added
    mutation_type_output=strings(ending_snp-starting_snp+1,1);
    nonsense_mutation_state=strings(ending_snp-starting_snp+1,1); %newly added
    essentiality=strings(ending_snp-starting_snp+1,1);
    non_evolvable_essentiality=strings(ending_snp-starting_snp+1,1);
    rm_af=repmat(-1,ending_snp-starting_snp+1,1);
    maf=repmat(-1,ending_snp-starting_snp+1,1);
    intergenic_gene_name_output=strings(ending_snp-starting_snp+1,1);
    chr_no_col=repmat(-1,ending_snp-starting_snp+1,1);
    chr_pos_col=repmat(-1,ending_snp-starting_snp+1,1);
    chr_pos_plotting_col=repmat(-1,ending_snp-starting_snp+1,1);

    for j=starting_snp+1:ending_snp+1 %+1 because starting snp and ending snp are 0 indexed, but the row numbers j of table X are 1-indexed
        j %comment out when running
        snp_index_number_array = split(X{j,2}{:,:},'snp');
        snp_index_number_array_size = size(snp_index_number_array);
        if snp_index_number_array_size(1,1)~=2
            fprintf('Chromosome end, no SNP-index at j=%d\n',j);
            gene_name_output(j-starting_snp,1)='END';
            gene_id_output(j-starting_snp,1)='END';
            mutation_type_output(j-starting_snp,1)='END';
            nonsense_mutation_state(j-starting_snp,1)=0; %newly added
            essentiality(j-starting_snp,1)='END';
            non_evolvable_essentiality(j-starting_snp,1)='END';
            rm_af(j-starting_snp,1)=-1;
            maf(j-starting_snp,1)=-1;
            intergenic_gene_name_output(j-starting_snp,1)='END';
            chr_no_col(j-starting_snp,1)=-1;
            chr_pos_col(j-starting_snp,1)=-1;
            chr_pos_plotting_col(j-starting_snp,1)=-1;
        end
        if snp_index_number_array_size(1,1)==2 %this skips all chromosome ends, where snp_index_number entry is 'END'
            snp_index_number = str2double(snp_index_number_array{2,1})-1; 
            
            %j
            [gene_name_output(j-starting_snp,1),gene_id_output(j-starting_snp,1)]=gene_name(snp_index_number); 


            if strcmp(gene_name_output(j-starting_snp, 1), 'N/A')
                mutation_type_output(j-starting_snp, 1) = 'intergenic';
                nonsense_mutation_state(j-starting_snp, 1) = '0';
            else
                gene_names_split = split(gene_name_output(j-starting_snp, 1), ';');
                gene_ids_split   = split(gene_id_output(j-starting_snp, 1), ';');
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
                mutation_type_output(j-starting_snp, 1) = mutation_type_temp;
                nonsense_mutation_state(j-starting_snp, 1) = nonsense_temp;
            end

            %checked both common and standard (YDR... format) gene names
            %against the list of essential and non-evolvable genes. The
            %table I used does not have a list of non-essential or NOT
            %non-evolvable genes. Can look for it later if needed, for now
            %non-essential and NOT non-evolvable can be annotated based on
            %just absensce in the above lists.

            %NEW BLOCK OF CODE, THAT CAN HANDLE MULTIPLE GENE NAMES MAPPING
            %TO A SINGLE SNP
            gene_names_split = split(gene_name_output(j-starting_snp, 1), ';');
            essentiality_temp = '';
            if strcmp(gene_name_output(j-starting_snp, 1), 'N/A')
                essentiality(j-starting_snp, 1) = 'N/A';
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
                  essentiality(j-starting_snp, 1) = essentiality_temp;
            end
            
            %NEW BLOCK OF CODE, THAT CAN HANDLE MULTIPLE GENE NAMES MAPPING
            %TO A SINGLE SNP
            gene_names_split = split(gene_name_output(j-starting_snp, 1), ';');
            non_evolvable_temp = '';
            if strcmp(gene_name_output(j-starting_snp, 1), 'N/A')
                non_evolvable_essentiality(j-starting_snp, 1) = 'N/A';
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
                   non_evolvable_essentiality(j-starting_snp, 1) = non_evolvable_temp;
            end
                
            snp_id=strcat('snp',num2str(snp_index_number+1));%+1 to go from 0 indexing to the 1 indexing used in the SNP number in BYxRM SNPs and all_allele_freqs
            l=find(strcmp(snp_id,X{:,2})); 
            corr_curr_QTL_idx=X{l,1}; %index of qtl in the 0-41609 range indexing used in all_allele_freqs.txt
            
            m=find(C{:,1}==corr_curr_QTL_idx);
            
            if isempty(m)~=1 %if we find a match in the 1011 genomes SNPs
                rm_af(j-starting_snp,1)=C{m,6};
                maf(j-starting_snp,1)=min(C{m,6},1-C{m,6}); %minor allele frequency
            end
            
            if isempty(m)==1 %if we do not find a match in the 1011 genomes SNPs
                rm_af(j-starting_snp,1)=-1; %some illogical value for allele frequency
                maf(j-starting_snp,1)=-1;
            end
    
            if strcmp(gene_name_output(j-starting_snp,1),'N/A')==1
                intergenic_gene_name_output(j-starting_snp,1)=intergenic_gene_name(snp_index_number);
                %j
                %intergenic_gene_name_output(j-starting_snp,1)
            else
                intergenic_gene_name_output(j-starting_snp,1)='genic';
            end
    

            %NEW BLOCK OF CODE, FIGURES OUT INTRONIC STATUS FOR ALL SNPS
            temp = intronic_status(snp_index_number);
            if ~strcmp(temp, 'N/A')
                intronic_gene = extractAfter(temp, 'intronic to ');
                
                current_gene_names = split(gene_name_output(j-starting_snp, 1), ';');
                is_in_gene_name = ismember(intronic_gene, current_gene_names);

                if is_in_gene_name
                    remove_idx = find(strcmp(current_gene_names, intronic_gene));

                    current_gene_names(remove_idx) = [];
                    if isempty(current_gene_names)
                        gene_name_output(j-starting_snp, 1) = 'N/A';
                    else
                        gene_name_output(j-starting_snp, 1) = strjoin(current_gene_names, ';');
                    end

                    current_gene_ids = split(gene_id_output(j-starting_snp, 1), ';');
                    current_gene_ids(remove_idx) = [];
                    if isempty(current_gene_ids)
                        gene_id_output(j-starting_snp, 1) = 'N/A';
                    else
                        gene_id_output(j-starting_snp, 1) = strjoin(current_gene_ids, ';');
                    end

                    current_essentiality = split(essentiality(j-starting_snp, 1), ';');
                    current_essentiality(remove_idx) = [];
                    if isempty(current_essentiality)
                        essentiality(j-starting_snp, 1) = 'N/A';
                    else
                        essentiality(j-starting_snp, 1) = strjoin(current_essentiality, ';');
                    end

                    current_non_evolvable = split(non_evolvable_essentiality(j-starting_snp, 1), ';');
                    current_non_evolvable(remove_idx) = [];
                    if isempty(current_non_evolvable)
                        non_evolvable_essentiality(j-starting_snp, 1) = 'N/A';
                    else
                        non_evolvable_essentiality(j-starting_snp, 1) = strjoin(current_non_evolvable, ';');
                    end

                    current_mutation_state = split(mutation_type_output(j-starting_snp, 1), ';');
                    current_mutation_state(remove_idx) = [];
                    if isempty(current_mutation_state)
                        mutation_type_output(j-starting_snp, 1) = 'intergenic';
                    else
                        mutation_type_output(j-starting_snp, 1) = strjoin(current_mutation_state, ';');
                    end

                    current_nonsense_state = split(nonsense_mutation_state(j-starting_snp, 1), ';');
                    current_nonsense_state(remove_idx) = [];
                    if isempty(current_nonsense_state)
                        nonsense_mutation_state(j-starting_snp, 1) = '0';
                    else
                        nonsense_mutation_state(j-starting_snp, 1) = strjoin(current_nonsense_state, ';');
                    end
                else
                    % case 4: unexpected
                    warning('intronic_status: gene %s from intronic annotation not found in gene_name_output (SNP index %d)', intronic_gene, snp_index_number);
                end  % end of if is_in_gene_name
                
                % append to intergenic_gene_name_output EXCEPT if SNP still belongs to an ORF - keep 'genic', do not append intronic annotation
                if strcmp(intergenic_gene_name_output(j-starting_snp, 1), 'genic') && ~strcmp(gene_name_output(j-starting_snp, 1), 'N/A')
                    % SNP still belongs to an ORF - keep 'genic', do not
                    % append intronic annotation %I think the first
                    % condition above is redundant, but keeping it
                elseif strcmp(intergenic_gene_name_output(j-starting_snp, 1), 'N/A') || strcmp(intergenic_gene_name_output(j-starting_snp, 1), 'genic')
                    intergenic_gene_name_output(j-starting_snp, 1) = temp;
                else
                    intergenic_gene_name_output(j-starting_snp, 1) = append(intergenic_gene_name_output(j-starting_snp, 1), ';', temp);
                end

                % append to intergenic_gene_name_output regardless of case
                %BELOW IS THE OLD BLOCK OF CODE. ONLY SMALL ISSUE: IT
                %ANNOTATES SNPS WITH THE INTERGENIC_GENE_NAME WITH
                %'INTRONIC TO ...' EVEN IF THEY ARE GENIC IN ANOTHER ORF
                %TO KEEP THINGS CONSISTENT WITH THE REST OF THE CODE, ALL
                %SNPS GENIC IN AN ORF HAVE THE INTERGENIC_GENE_NAME
                %'GENIC', THATS IT
                %ABOVE IS A NEW  BLOCK OF CODE. NEEDS TO BE RUN ON ALL SNPS
                %AGAIN IF NEEDED, I MANUALLY EDITED THE DATA FOR THE 4
                %TOTAL SNPS THAT ARE IN ONE ORF BUT INTRONIC TO ANOTHER
                % if strcmp(intergenic_gene_name_output(j-starting_snp, 1), 'N/A') || strcmp(intergenic_gene_name_output(j-starting_snp, 1), 'genic') %or genic in case genic status is removed because it is intronic
                %     intergenic_gene_name_output(j-starting_snp, 1) = temp;
                % else
                %     intergenic_gene_name_output(j-starting_snp, 1) = append(intergenic_gene_name_output(j-starting_snp, 1), ',', temp);
                % end
            end  % end of if ~strcmp(temp, 'N/A')
            %END OF INTRONIC STATUS BLOCK


            if strcmp(gene_name_output(j-starting_snp,1),'N/A')==1 & strcmp(mutation_type_output(j-starting_snp,1),'intergenic')==0
                fprintf('no gene name is annotated but SNP is also not annotated as intergenic at SNP index = %d\n',j);
            end
            
            [chr_no_col(j-starting_snp,1),chr_pos_col(j-starting_snp,1)]=gene_pos_chr_idx(snp_index_number);
            chr_pos_plotting_col(j-starting_snp,1)=chr_pos_plotting(chr_no_col(j-starting_snp,1),chr_pos_col(j-starting_snp,1));
            %j
        end
    end
    X_subset = X(starting_snp+1:ending_snp+1,:);
    
    X_subset.("gene_name")=gene_name_output;
    X_subset.("mutation_type")=mutation_type_output;
    X_subset.("gene_essentiality")=essentiality;
    X_subset.("non_evolvable_gene_essentiality")=non_evolvable_essentiality;
    X_subset.("minor_allele_frequency")=maf;
    X_subset.("rm_allele_frequency")=rm_af;
    X_subset.("intergenic_gene_name")=intergenic_gene_name_output;
    X_subset.("chr_no")=chr_no_col;
    X_subset.("chr_pos")=chr_pos_col;
    X_subset.("chr_pos_plotting")=chr_pos_plotting_col;
    X_subset.("nonsense_mutation_state")=nonsense_mutation_state; %newly added
    X_subset.("gene_id")=gene_id_output;

    table_name=strcat('output/BYxRM_SNPs_detailed_annotated_part_',num2str(i),'.txt');
    writetable(X_subset,table_name,Delimiter='\t');
    warning(warnState);
end
