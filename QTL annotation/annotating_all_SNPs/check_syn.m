function [name1, name2] = check_syn(gene_pos_bp,gene_pos_chr,gene_start,gene_stop,by_allele,rm_allele,gene_strand,gene_pos)
%quick mention of strand information in gff files: encoded as + and -
%top strand genes can only be transcribed "forwards" (left to right) and
%bottom strand genes can only be transcribed "reverse" (right to left), so
%one variable is sufficient to encode both top vs bottom strand and forward
%vs reverse orientation of genes. This is because RNA can only be made in
%the 5' -> 3' direction

ref_genome=fastaread('RM11a_bbq_github.fasta'); 
ref_chr=ref_genome(gene_pos_chr).Sequence; 
%chromosome of interest (needs to be gene_pos_chr+1 for the BY ref genome. 
% Note that the by ref genome has a large ~20kb truncation at the end of chr13, 
% prob an assembly error
%based on length, chr17 is mitochondria for the BY ref genome

%the RM ref genome doesn't have a 2-micron chr assembly or a mito genome
%assembly. Not an issue really because none of the 42k variants are called
%to be mito or 2-micron genomic.


%need to switch from ascending to descending order for - strand, multi-exon
%genes
if strcmp(gene_strand, '-')
    gene_start = fliplr(gene_start);
    gene_stop  = fliplr(gene_stop);
end

% concatenate all CDS fragments in order
ref_region = '';
for k = 1:length(gene_start)
    ref_region = [ref_region, ref_chr(gene_start(k):gene_stop(k))];
end

% find which CDS fragment contains gene_pos_bp and what offset it is at
cumulative_offset = 0;
snp_pos_in_concat = NaN;
for k = 1:length(gene_start)
    if gene_pos_bp >= gene_start(k) && gene_pos_bp <= gene_stop(k)
        snp_pos_in_concat = cumulative_offset + (gene_pos_bp - gene_start(k) + 1);
        break
    end
    cumulative_offset = cumulative_offset + (gene_stop(k) - gene_start(k) + 1);
end

if isnan(snp_pos_in_concat)
    error('check_syn: gene_pos_bp %d not found in any CDS fragment', gene_pos_bp)
end

by_region = ref_region;
rm_region = ref_region;
by_region(snp_pos_in_concat) = by_allele{:,:};
rm_region(snp_pos_in_concat) = rm_allele{:,:};

% fprintf('DEBUG: gene_strand is "%s", length %d\n', gene_strand, length(gene_strand))
% fprintf('DEBUG: ref_region length = %d\n', length(ref_region))
% fprintf('DEBUG: first 10 nt of ref_region = %s\n', ref_region(1:10))
% 
% fprintf('DEBUG: strcmp result for - : %d\n', strcmp('-', gene_strand))
% fprintf('DEBUG: strcmp result for + : %d\n', strcmp('+', gene_strand))


%note than ref protein here will be the RM coding sequence as the genome 
%reference fasta used is the RM version
if strcmp('-',gene_strand)==1
    by_protein=nt2aa(seqrcomplement(by_region));
    rm_protein=nt2aa(seqrcomplement(rm_region));
    ref_protein=nt2aa(seqrcomplement(ref_region));
elseif strcmp('+',gene_strand)==1
    by_protein=nt2aa(by_region);
    rm_protein=nt2aa(rm_region);
    ref_protein=nt2aa(ref_region);
end

%fprintf('DEBUG: rm_protein seq is = %s\n', rm_protein)

% ── Check ref_protein starts with M (start codon, methionine) and ends with * (stop codon) ─────────────────────────
if ref_protein(1) ~= 'M' || ref_protein(end) ~= '*'
    fprintf('ERROR: ref_protein does not start with Met or end with a stop codon.\n')
    fprintf('  gene_pos, 0-indexed:  %d\n', gene_pos)
    fprintf('  gene_pos_bp:  %d\n', gene_pos_bp)
    fprintf('  gene_pos_chr: %d\n', gene_pos_chr)
    fprintf('  gene_start:   %d\n', gene_start)
    fprintf('  gene_stop:    %d\n', gene_stop)
    fprintf('  by_allele:    %s\n', by_allele{:,:})
    fprintf('  rm_allele:    %s\n', rm_allele{:,:})
    fprintf('  gene_strand:  %s\n', gene_strand)
    fprintf('  ref_protein:  %s\n', ref_protein)
    warning('check_syn: invalid ref_protein — does not begin with M (%s) or end with * (%s)', ...
          ref_protein(1), ref_protein(end))
end

if strcmp(rm_protein,ref_protein)==0 %rm_protein does not perfectly match ref_protein
    fprintf('ERROR: ref_protein does not match rm_protein sequence.\n')
    fprintf('  ref_protein:  %s\n', ref_protein)
    fprintf('  rm_protein:  %s\n', rm_protein)
    warning('check_syn: ref_protein and rm_protein sequences do not match')
end

% initialize outputs
name1 = '';
name2 = 0;

if strcmp(by_protein, rm_protein) == 0
    name1 = 'non-synonymous';
    diff_pos = find(by_protein ~= rm_protein, 1, 'first');
    if ~isempty(diff_pos) && by_protein(diff_pos) == '*'
        name2 = 1;
    end
elseif strcmp(by_protein, rm_protein) == 1
    name1 = 'synonymous';
end
end
