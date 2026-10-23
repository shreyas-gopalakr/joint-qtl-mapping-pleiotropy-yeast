A=readtable('final_model_cleaned_linear_effects_annotated_v8.txt');
C=readtable('pos_uncerts.csv');

A.Var1(2:height(A))==C{:,1}; %all 1s, so positions match across both arrays

C_row1 = table(-1, -1, -1, 'VariableNames', {'lead', 'bottom', 'top'});
C = [C_row1; C];

A.CI_lead_SNP=C{:,1};
A.CI_bottom_SNP=C{:,2};
A.CI_top_SNP=C{:,3};

A.CI_genes(1)="CI_genes";
for i=2:height(A)
    A.CI_genes(i)=CI_genes(A.CI_bottom_SNP(i),A.CI_top_SNP(i));
end

writetable(A,'final_model_cleaned_linear_effects_annotated_v9_temp.txt',Delimiter='\t');

nnz(A.CI_genes=="N/A"); %4 QTLs whose CIs just say N/A for now
%To do: 1. modify the intergenic_gene_name.m code to check if any genes are
%missing whose intergenic regions are not included in the current list
%2. Also run the intronic gene code and check things

CI_size_distribution=count(A.CI_genes,',')+1;
histogram(CI_size_distribution,BinEdges=[1:max(CI_size_distribution)]);
%only 35 QTLs mapped to single-gene precision