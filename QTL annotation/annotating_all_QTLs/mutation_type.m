function [name1, name2] = mutation_type(gene_pos, idx, gene_name, gene_id)
%modified to annotate the mutation_types for all genes/CDS' that the SNP
%can be associated with
X = readtable('BYxRM_SNPs_detailed.xlsx');

if idx == 0
    idx2 = 1;
elseif idx == 1
    idx2 = 0;
end

match_row = find(strcmp(X{:,2}, strcat('snp', num2str(gene_pos + idx2))), 1, 'first');
if isempty(match_row)
    name1 = 'error: no chromosomal position matches this SNP index';
    name2 = 'error: no chromosomal position matches this SNP index';
    return
end
gene_pos_bp  = X{match_row, 4};
gene_pos_chr = X{match_row, 3};
by_allele    = X{match_row, 5};
rm_allele    = X{match_row, 6};

V = GFFAnnotation('BY4742_sorted_annotations.gff3');

for j = 1:V.NumEntries
    gene_feature = V.Feature(j);
    ref_chr_temp = V.Reference(j);
    ref_chr_temp2 = ref_chr_temp{1,1};
    ref_chr = str2num(ref_chr_temp2(end-1:end));

    if gene_pos_bp >= V.Start(j) && gene_pos_bp <= V.Stop(j)
        if strcmp(gene_feature{1,1}, 'CDS') == 1
            if ref_chr == gene_pos_chr
                pid = parse_parent(V.Attributes{j});
                if strcmp(pid, gene_id) || is_mrna_child(V, pid, gene_id)
                    gene_start  = V.Start(j);
                    gene_stop   = V.Stop(j);
                    gene_strand = V.Strand(j);
                    gene_strand = gene_strand{:};
                    parent_id   = pid;
                    name1       = 'genic';
                    break
                end
            end
        end
    end
end

% ── 1. Check gene_strand is valid ────────────────────────────────────────────
if exist("name1") == 1 && strcmp(name1, 'genic')
    if ~strcmp(gene_strand, '+') && ~strcmp(gene_strand, '-')
        fprintf('ERROR: gene_strand is neither + nor -\n')
        fprintf('  gene_pos:     %d\n', gene_pos)
        fprintf('  gene_pos_bp:  %d\n', gene_pos_bp)
        fprintf('  gene_pos_chr: %d\n', gene_pos_chr)
        fprintf('  gene_start:   %d\n', gene_start)
        fprintf('  gene_stop:    %d\n', gene_stop)
        fprintf('  gene_strand:  %s\n', gene_strand)
        error('mutation_type: gene_strand is invalid: "%s"', gene_strand)
    end
end

% ── 2. Search for multiple CDS / intron features with the same parent ─────────
if exist("name1") == 1 && strcmp(name1, 'genic')

    cds_starts  = [];
    cds_stops   = [];
    intron_found = false;
    frameshift_found = false;

    for j = 1:V.NumEntries
        feat = V.Feature(j);
        par  = parse_parent(V.Attributes{j});

        if ~strcmp(par, parent_id)
            continue  % only look at features sharing the same parent mRNA
        end

        if strcmp(feat{1,1}, 'CDS')
            cds_starts(end+1) = V.Start(j);
            cds_stops(end+1)  = V.Stop(j);
        elseif strcmp(feat{1,1}, 'intron')
            intron_found = true;
        elseif strcmp(feat{1,1}, 'plus_1_translational_frameshift')
            frameshift_found = true;
        end
    end

    % sort CDS fragments by coordinate
    [cds_starts, sort_idx] = sort(cds_starts);
    cds_stops = cds_stops(sort_idx);

    num_cds = length(cds_starts);

    if num_cds > 1 && ~intron_found && ~frameshift_found
        fprintf('ERROR: multiple CDS features found for parent %s but no intron or frameshift feature.\n', parent_id)
        fprintf('  gene_pos:     %d\n', gene_pos)
        fprintf('  gene_pos_bp:  %d\n', gene_pos_bp)
        fprintf('  gene_pos_chr: %d\n', gene_pos_chr)
        fprintf('  gene_strand:  %s\n', gene_strand)
        fprintf('  CDS starts:   %s\n', num2str(cds_starts))
        fprintf('  CDS stops:    %s\n', num2str(cds_stops))
        error('mutation_type: multiple CDS features but no intron or frameshift for parent %s', parent_id)
    end

    % For frameshift genes, translate only the CDS fragment containing the SNP
    % (each CDS fragment is in a different reading frame, so they must be
    % translated independently rather than concatenated)
    if num_cds > 1 && frameshift_found
        snp_cds_idx = NaN;
        for k = 1:num_cds
            if gene_pos_bp >= cds_starts(k) && gene_pos_bp <= cds_stops(k)
                snp_cds_idx = k;
                break
            end
        end
        if isnan(snp_cds_idx)
            error('mutation_type: SNP at %d not in any CDS fragment for frameshift gene %s', gene_pos_bp, parent_id)
        end
        cds_starts = cds_starts(snp_cds_idx);
        cds_stops  = cds_stops(snp_cds_idx);
    end

    % pass arrays of starts/stops to check_syn (single CDS is just a 1-element array)
    %HAVING ISSUES IN CERTAIN PARTS OF CERTAIN CHROMOSOMES (EG: END OF CHR
    %13)
    %BEST TO TRANSLATE ALL ORFS FROM BY4742.GFF ANNOTATIONS TO THE STANDARD
    %S288C-R64-1-1.GFF ANNOTATIONS AND THEN USE THE S288C-R64-1-1.FASTA FOR
    %THE SEQUENCE OF THE PROTEINS
    %CAN DOUBLE CHECK THE TYPE OF MUTATION CALLED BY EACH

    %DON'T NEED TO DO THE ABOVE, FOUND THE REF FASTAS
    %BY REF IS MISSING THE END OF CHR13, RM SEEMS COMPLETE
    %USE THE RM REF, CHECK THAT THE REF INDEED CONTAINS THE CORRECT ALLELE
    %AT THE ADVERTISED LOCATION
    [name1, name2] = check_syn(gene_pos_bp, gene_pos_chr, cds_starts, cds_stops, by_allele, rm_allele, gene_strand, gene_pos);

elseif exist("name1") == 0
    name1 = 'intergenic';
    name2 = 0;
end

end

function parent = parse_parent(attributes_str)
% attributes_str looks like:
% 'Parent=YDR025W_mRNA;Name=YDR025W_CDS;orf_classification=Verified'
parent = '';
parts = strsplit(attributes_str, ';');   % split on ; to get individual key=value pairs
for k = 1:length(parts)
    kv = strsplit(parts{k}, '=');        % split each pair on =
    if length(kv) == 2 && strcmp(strtrim(kv{1}), 'Parent')
        parent = strtrim(kv{2});
        return
    end
end
end

function result = is_mrna_child(V, mrna_id, gene_id)
% checks if mrna_id is an mRNA whose Parent is gene_id
result = false;
for j = 1:V.NumEntries
    feat = V.Feature(j);
    if strcmp(feat{1,1}, 'mRNA')
        % check if this mRNA's ID matches mrna_id
        attr_str = V.Attributes{j};
        parts = strsplit(attr_str, ';');
        this_id = '';
        this_parent = '';
        for k = 1:length(parts)
            kv = strsplit(parts{k}, '=');
            if length(kv) == 2
                if strcmp(strtrim(kv{1}), 'ID')
                    this_id = strtrim(kv{2});
                elseif strcmp(strtrim(kv{1}), 'Parent')
                    this_parent = strtrim(kv{2});
                end
            end
        end
        if strcmp(this_id, mrna_id) && strcmp(this_parent, gene_id)
            result = true;
            return
        end
    end
end
end