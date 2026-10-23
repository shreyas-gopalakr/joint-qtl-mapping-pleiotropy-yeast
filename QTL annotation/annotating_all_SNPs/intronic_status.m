function name1 = intronic_status(gene_pos)
%all the following text is copied over from the function gene_name.m
    X=readtable('BYxRM_nanopore_SNPs.txt','Delimiter','\t'); %list of SNP positions
    Y=table2array(X); %convert table to array
    Z=Y(Y~=2); %removing 2's, ends of chr are stored as 2 (stupid)
    Z_size=size(Z);
    chr(1,1)=1;
    count=2; %count of number of chrs
    for i=1:Z_size(2)-1
        if Z(i)>Z(i+1) %kind of a stupid way of finding out where a chr starts, but works here. Can just use position of the 2's in future.
            chr(1,count)=i+1; %chr stores the positions of the start of a chromosome in SNP #
            count=count+1;
        end
    end
    chr(1,count)=Z_size(2)+1; %end of copying from lod_score_by_env_chr_number.m script
    gene_pos_bp=Z(gene_pos+1); %input data is 0 indexed, hence the +1 here
    
    chr_num=size(chr);
    for i=1:chr_num(1,2)-1
        if gene_pos+1>=chr(1,i) & gene_pos+1<chr(1,i+1)
            gene_pos_chr=i;
        end
    end
    
    
    V=GFFAnnotation('BY4742_sorted_annotations.gff3');


    for j=1:V.NumEntries
    gene_feature=V.Feature(j);
    ref_chr_temp=V.Reference(j);
    ref_chr_temp2=ref_chr_temp{1,1};
    ref_chr=str2num(ref_chr_temp2(end-1:end));
    if gene_pos_bp>=V.Start(j) & gene_pos_bp<=V.Stop(j)==1 %within the feature at j
        if strcmp(gene_feature{1,1},"intron")==1 %the feature is an intron
            if ref_chr==gene_pos_chr %feature is in the correct chromosome
                % check parent gene orf_classification - skip dubious ORFs
                attr = V.Attributes(j);
                attr_str = attr{1,1};
                attr_split = split(attr_str, ';');
                
                [orf_class, gene_name] = get_parent_gene_orf_class(V, attr_str);
                if strcmp(orf_class, 'Dubious')
                    continue
                end
                
                if isempty(gene_name)
                    warning('intronic_status: could not find gene name for intron at gene_pos %d', gene_pos);
                    gene_name = 'unknown';
                end
                name1 = append('intronic to ', gene_name);

            end
        end
    end
    end

    if exist('name1')~=1 %if no name is assigned, because it is intergenic or whatever
    name1='N/A';
    end
end


function [orf_class, gene_name] = get_parent_gene_orf_class(V, intron_attr_str)
orf_class = '';
gene_name = '';
% Given an intron's attribute string, traverses up to the parent gene
% and returns its orf_classification. Returns '' if not found.

% step 1: get the intron's Parent (e.g. YBR048W_mRNA)
parts = strsplit(intron_attr_str, ';');
mrna_id = '';
for k = 1:length(parts)
    kv = strsplit(parts{k}, '=');
    if length(kv) == 2 && strcmp(strtrim(kv{1}), 'Parent')
        mrna_id = strtrim(kv{2});
        break
    end
end
if isempty(mrna_id)
    return
end

% step 2: find intermediate feature with ID=mrna_id
gene_id = '';
nc_types = {'mRNA', 'tRNA', 'snRNA', 'rRNA', 'snoRNA', 'ncRNA'};
for j = 1:V.NumEntries
    feat = V.Feature(j);
    if ismember(feat{1,1}, nc_types)
        attr_str = V.Attributes{j};
        parts2 = strsplit(attr_str, ';');
        this_id = '';
        for k = 1:length(parts2)
            kv = strsplit(parts2{k}, '=');
            if length(kv) == 2 && strcmp(strtrim(kv{1}), 'ID')
                this_id = strtrim(kv{2});
                break
            end
        end
        if strcmp(this_id, mrna_id)
            if strcmp(feat{1,1}, 'mRNA')
                % normal case: get Parent gene_id and continue to step 3
                for k = 1:length(parts2)
                    kv = strsplit(parts2{k}, '=');
                    if length(kv) == 2 && strcmp(strtrim(kv{1}), 'Parent')
                        gene_id = strtrim(kv{2});
                        break
                    end
                end
            else
                % tRNA/snRNA etc: this IS the top-level feature, extract name directly
                orf_class = 'ncRNA';
                for k = 1:length(parts2)
                    kv = strsplit(parts2{k}, '=');
                    if length(kv) == 2 && strcmp(strtrim(kv{1}), 'Name')
                        gene_name = strtrim(kv{2});
                        break
                    end
                end
                return  % done, no step 3 needed
            end
            break
        end
    end
end

% step 3: find the gene feature with ID=gene_id and get its orf_classification
for j = 1:V.NumEntries
    feat = V.Feature(j);
    if strcmp(feat{1,1}, 'gene')
        attr_str = V.Attributes{j};
        parts3 = strsplit(attr_str, ';');
        this_id = '';
        for k = 1:length(parts3)
            kv = strsplit(parts3{k}, '=');
            if length(kv) == 2 && strcmp(strtrim(kv{1}), 'ID')
                this_id = strtrim(kv{2});
                break
            end
        end
        if strcmp(this_id, gene_id)
            % found the gene, now get orf_classification and gene name
            for k = 1:length(parts3)
                if contains(parts3{k}, 'orf_classification=')
                    orf_temp = split(parts3{k}, 'orf_classification=');
                    orf_class = orf_temp{2,1};
                end
                if contains(parts3{k}, 'gene=')
                    gene_temp = split(parts3{k}, 'gene=');
                    gene_name = gene_temp{2,1};
                end
            end
            if isempty(gene_name)
                % fall back to Name= field on the gene feature
                for k = 1:length(parts3)
                    if contains(parts3{k}, 'Name=')
                        name_temp = split(parts3{k}, 'Name=');
                        gene_name = name_temp{2,1};
                        break
                    end
                end
            end
            return
        end
    end
end
end