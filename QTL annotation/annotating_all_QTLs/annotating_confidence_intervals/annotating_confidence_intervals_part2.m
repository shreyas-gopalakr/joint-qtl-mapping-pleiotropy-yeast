A=readtable('final_model_cleaned_linear_effects_annotated_v9_temp.txt');

A.CI_genes_with_intergenic(1)="CI_genes_with_intergenic";
for i=2:height(A)
    i
    A.CI_genes_with_intergenic(i)=CI_genes_with_intergenic(A.CI_bottom_SNP(i),A.CI_top_SNP(i),A.CI_genes(i));
end

writetable(A,'final_model_cleaned_linear_effects_annotated_v9.txt',Delimiter='\t');

