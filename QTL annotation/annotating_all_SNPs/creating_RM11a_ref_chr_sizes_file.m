A=fastaread("RM11a_bbq_github.fasta");
B=readtable('Scer_genome_ref_chr_sizes.tsv',"FileType","text",'Delimiter', '\t');

len_chr=zeros(height(A),1);
for i=1:height(A)
    chr=A(i).Sequence;
    len_chr(i,1)=length(chr);
end

V = GFFAnnotation('BY4742_sorted_annotations.gff3');
idx = strcmp(V.Feature, 'contig');

% extract fields manually using the logical index
contig_names = V.Reference(idx);
contig_starts = V.Start(idx);
contig_stops = V.Stop(idx);

%contig starts are all 1, so the length of the chr/contig is just the 
%value of contig_stops

chr_lengths_gff = contig_stops(2:17,1); %removing row 1 and 18, correspond
%to 2-micron and mito, which doesn't exist in the RM11a ref fasta


chr_lengths_equality=(chr_lengths_gff==len_chr);
%all chr lengths are equal from both the gff and associated fasta, except
%chr13, where the lengths differ by 1bp

%chr_lengths_gff has the 1bp larger size for chr13, so I'll just use that
%for now

contig_names_trimmed = contig_names(2:17, 1);

T = table(contig_names_trimmed, B.ChromosomeName(1:16), B.UCSCStyleName(1:16), chr_lengths_gff, 'VariableNames', {'contig_name', 'ChromosomeName', 'UCSCStyleName', 'SeqLength'});

writetable(T, 'RM11a_ref_chr_sizes.tsv', 'FileType', 'text', 'Delimiter', '\t');
