%largely copied over from a function called gene_name.m from /Users/shreyasg/Downloads/EF_joint_inference_Jul_2025/2 - Gene name annotation
function CI_genes_output = CI_genes(CI_bottom_SNP,CI_top_SNP)
CI_genes_output="";
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
CI_bottom_bp=Z(CI_bottom_SNP+1); %input data is 0 indexed, hence the +1 here
CI_top_bp=Z(CI_top_SNP+1);

%finding chr number for CI bottom SNP
chr_num=size(chr);
for i=1:chr_num(1,2)-1
    if CI_bottom_SNP+1>=chr(1,i) & CI_bottom_SNP+1<chr(1,i+1)
        CI_bottom_chr=i;
    end
end

%finding chr number for CI top SNP
for i=1:chr_num(1,2)-1
    if CI_top_SNP+1>=chr(1,i) & CI_top_SNP+1<chr(1,i+1)
        CI_top_chr=i;
    end
end

V=GFFAnnotation('BY4742_sorted_annotations.gff3');

if CI_top_chr~=CI_bottom_chr
    warning('The top and bottom SNPs correspond to different chromosomes')
end

for j=1:V.NumEntries
    gene_feature=V.Feature(j);
    ref_chr_temp=V.Reference(j);
    ref_chr_temp2=ref_chr_temp{1,1};
    ref_chr=str2num(ref_chr_temp2(end-1:end));
    cond1=CI_bottom_bp>=V.Start(j) & CI_bottom_bp<=V.Stop(j); %gene containing bottom SNP
    cond2=CI_top_bp>=V.Start(j) & CI_top_bp<=V.Stop(j); %gene containing top SNP
    cond3=CI_bottom_bp<V.Start(j) & CI_bottom_bp<V.Stop(j) & CI_top_bp>V.Start(j) & CI_top_bp>V.Stop(j); %gene between bottom and top
    if (cond1 | cond2 | cond3) %satisfying one of the 3 conditions mentioned above
        if strcmp(gene_feature{1,1},"gene")==1 %the feature is a gene
            if ref_chr==CI_top_chr %feature is in the correct chromosome
                attr=V.Attributes(j);
                attr_str=attr{1,1};
                attr_split=split(attr_str,';');
                attr_size=size(attr_split);

                % check orf_classification - skip dubious ORFs
                orf_class = '';
                for k = 1:length(attr_split)
                    if contains(attr_split{k,1}, 'orf_classification=')
                        orf_temp = split(attr_split{k,1}, 'orf_classification=');
                        orf_class = orf_temp{2,1};
                        break
                    end
                end
                if strcmp(orf_class, 'Dubious')
                    continue  % skip this gene feature entirely
                end

                    if (CI_genes_output == "") %adding the first gene of CI
                        if contains(attr_split{4,1},'gene=')==1
                            name1_temp=split(attr_split{4,1},'gene=');
                            CI_genes_output=name1_temp{2,1};
                        else 
                            name1_temp=split(attr_split{2,1},'Name=');
                            CI_genes_output=name1_temp{2,1};
                        end
                    else %adding the second or later gene of CI
                        if contains(attr_split{4,1},'gene=')==1
                            name1_temp=split(attr_split{4,1},'gene=');
                            CI_genes_output=strcat(CI_genes_output,";",name1_temp{2,1});
                        else 
                            name1_temp=split(attr_split{2,1},'Name=');
                            CI_genes_output=strcat(CI_genes_output,";",name1_temp{2,1});                        
                        end
                    end
                %end
            end
        end
    end
end

if (CI_genes_output == "") %if no name is assigned, because the full region is intergenic or whatever
CI_genes_output='N/A';
end

end