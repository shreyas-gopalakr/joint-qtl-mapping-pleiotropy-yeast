function [name1, name2] = gene_name(gene_pos) %going to return gene name and gene "ID" (YDR... name), its useful in figuring out associated mRNA and CDS fields
%name1 returns gene names and name2 returns gene IDs
%edited to return multiple gene names if the SNP belongs in multiple genes
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

name1 = '';
name2 = ''; % will store comma-separated list of gene IDs


for j=1:V.NumEntries
    gene_feature=V.Feature(j);
    ref_chr_temp=V.Reference(j);
    ref_chr_temp2=ref_chr_temp{1,1};
    ref_chr=str2num(ref_chr_temp2(end-1:end));
    if gene_pos_bp>=V.Start(j) & gene_pos_bp<=V.Stop(j)==1 %within the feature at j
        if strcmp(gene_feature{1,1},"gene")==1 %the feature is a gene
            if ref_chr==gene_pos_chr %feature is in the correct chromosome
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


                if contains(attr_split{4,1}, 'gene=') == 1
                    name1_temp = split(attr_split{4,1}, 'gene=');
                    name1_temp2 = name1_temp{2,1};
                else
                    name1_temp = split(attr_split{2,1}, 'Name=');
                    name1_temp2 = name1_temp{2,1};
                end
                
                % extract ID field
                id_temp = split(attr_split{1,1}, 'ID=');
                id_temp2 = id_temp{2,1};
                
                if isempty(name1) %if this is the first matching gene name and id found
                    name1 = name1_temp2;
                    name2 = id_temp2;
                else %if multiple gene names and ids match this region
                    name1 = append(name1, ';', name1_temp2);
                    name2 = append(name2, ';', id_temp2);
                end

            end
        end
    end
end

%if no name is assigned, because it is intergenic or whatever
if isempty(name1)
    name1 = 'N/A';
    name2 = 'N/A';
end

end