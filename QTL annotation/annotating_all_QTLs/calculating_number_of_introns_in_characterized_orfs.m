%incomplete annotation at the intron level:
% Dubious: 3
% Uncharacterized: 37
% Verified: 271
% Verified|silenced_gene: 2
% not found: 63

%all the not found orf annotation cases correspond to tRNAs (>90%), snRNAs,
%etc

%63 introns don't have an orf classification status mapped
%Here was the result for gene features:
% Dubious: 785
% Uncharacterized: 839
% Verified: 4970
% Verified|silenced_gene: 4

V = GFFAnnotation('BY4742_sorted_annotations.gff3');
orf_classifications = {};
for j = 1:V.NumEntries
    gene_feature = V.Feature(j);
    if strcmp(gene_feature{1,1}, 'intron') == 1
        attr = V.Attributes(j);
        attr_str = attr{1,1};
        
        [orf_class, gene_name] = get_parent_gene_orf_class(V, attr_str);

        if isempty(orf_class)
            % print for debugging - only for the 63 failing cases
            attr_split = split(attr_str, ';');
            intron_name = 'unknown';
            for k = 1:length(attr_split)
                if contains(attr_split{k}, 'Name=')
                    name_temp = split(attr_split{k}, 'Name=');
                    intron_name = name_temp{2,1};
                    break
                end
            end
            fprintf('intron: %s | parent gene found: %s | orf_class: %s\n', intron_name, gene_name, orf_class);
            orf_class = 'not found even after parent traversal';
        end
        
        orf_classifications{end+1} = orf_class;
    end
end

% count occurrences of each unique classification
[unique_classes, ~, ic] = unique(orf_classifications);
for k = 1:length(unique_classes)
    fprintf('%s: %d\n', unique_classes{k}, sum(ic == k));
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

% step 2: find the mRNA feature with ID=mrna_id and get its Parent (gene ID)
gene_id = '';
for j = 1:V.NumEntries
    feat = V.Feature(j);
    if strcmp(feat{1,1}, 'mRNA')
        attr_str = V.Attributes{j};
        parts2 = strsplit(attr_str, ';');
        this_id = '';
        this_parent = '';
        for k = 1:length(parts2)
            kv = strsplit(parts2{k}, '=');
            if length(kv) == 2
                if strcmp(strtrim(kv{1}), 'ID')
                    this_id = strtrim(kv{2});
                elseif strcmp(strtrim(kv{1}), 'Parent')
                    this_parent = strtrim(kv{2});
                end
            end
        end
        if strcmp(this_id, mrna_id)
            gene_id = this_parent;
            break
        end
    end
end
if isempty(gene_id)
    return
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