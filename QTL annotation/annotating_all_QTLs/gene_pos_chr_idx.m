function [gene_pos_chr,gene_pos_bp] = gene_pos_chr_idx(gene_pos)
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
end