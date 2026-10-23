function chr_pos_plot_out = chr_pos_plotting(gene_pos_chr,gene_pos_bp)
%make a new table here based on our assembly's chr sizes
X=readtable('RM11a_ref_chr_sizes.tsv',"FileType","text",'Delimiter', '\t');
%sum(X.SeqLength(1:16)) %ans is 12071326
i=1;
chr_pos_plot_out=0;
    while i<gene_pos_chr
        chr_pos_plot_out=chr_pos_plot_out + X.SeqLength(i);
        i=i+1;
    end
chr_pos_plot_out=chr_pos_plot_out + gene_pos_bp;
end