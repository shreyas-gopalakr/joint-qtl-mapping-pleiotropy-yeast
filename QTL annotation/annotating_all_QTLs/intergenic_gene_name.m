function name1 = intergenic_gene_name(gene_pos)
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

for j=1:V.NumEntries
    gene_feature=V.Feature(j);
    ref_chr_temp=V.Reference(j);
    ref_chr_temp2=ref_chr_temp{1,1};
    ref_chr=str2num(ref_chr_temp2(end-1:end));
    %looking for upstream regulatory regions
    if strcmp(gene_feature{1,1}, "gene") == 1
        if ref_chr == gene_pos_chr
            gene_strand_temp = V.Strand(j);
            gene_strand_temp = gene_strand_temp{1,1};
        

        attr_temp = V.Attributes(j);
        attr_str_temp = attr_temp{1,1};
        attr_split_temp = split(attr_str_temp, ';');    
        %following step skips dubious orfs
        if any(contains(attr_split_temp, 'orf_classification=Dubious'))
            continue
        end

        if strcmp(gene_strand_temp, '+')
            is_upstream = gene_pos_bp < V.Start(j) && gene_pos_bp >= V.Start(j) - 500;
        elseif strcmp(gene_strand_temp, '-')
            is_upstream = gene_pos_bp > V.Stop(j) && gene_pos_bp <= V.Stop(j) + 500;
        else
            is_upstream = false;
        end
        
        if is_upstream
                attr=V.Attributes(j);
                attr_str=attr{1,1};
                attr_split=split(attr_str,';');
                attr_size=size(attr_split);
                if contains(attr_split{4,1},'gene=')==1
                    name1_temp=split(attr_split{4,1},'gene=');
                    name1_temp2=name1_temp{2,1};
                    if exist('name1','var')==0
                        name1=append("-",name1_temp2);
                    else
                        name1=append(name1,";-",name1_temp2);
                        %j
                    end
                else 
                    name1_temp=split(attr_split{2,1},'Name=');
                    name1_temp2=name1_temp{2,1};
                    if exist('name1','var')==0
                        name1=append("-",name1_temp2);
                    else
                        name1=append(name1,";-",name1_temp2);
                    end
                end
                end
            end
        end
     %looking for downstream regulatory regions
    if strcmp(gene_feature{1,1}, "gene") == 1
        if ref_chr == gene_pos_chr
            gene_strand_temp = V.Strand(j);
            gene_strand_temp = gene_strand_temp{1,1};
        

        attr_temp = V.Attributes(j);
        attr_str_temp = attr_temp{1,1};
        attr_split_temp = split(attr_str_temp, ';');    
        %following step skips dubious orfs
        if any(contains(attr_split_temp, 'orf_classification=Dubious'))
            continue
        end

        if strcmp(gene_strand_temp, '+')
            is_downstream = gene_pos_bp > V.Stop(j) && gene_pos_bp <= V.Stop(j) + 500;
        elseif strcmp(gene_strand_temp, '-')
            is_downstream = gene_pos_bp < V.Start(j) && gene_pos_bp >= V.Start(j) - 500;
        else
            is_downstream = false;
        end
        
        if is_downstream
                attr=V.Attributes(j);
                attr_str=attr{1,1};
                attr_split=split(attr_str,';');
                attr_size=size(attr_split);
                if contains(attr_split{4,1},'gene=')==1
                    name1_temp=split(attr_split{4,1},'gene=');
                    name1_temp2=name1_temp{2,1};
                    if exist('name1','var')==0
                        name1=append("+",name1_temp2);
                    else
                        name1=append(name1,";+",name1_temp2);
                    end
                else 
                    name1_temp=split(attr_split{2,1},'Name=');
                    name1_temp2=name1_temp{2,1};
                    if exist('name1','var')==0
                        name1=append("+",name1_temp2);
                    else
                        name1=append(name1,";+",name1_temp2);
                    end
                end
         end
        end
     end
end

if exist('name1')~=1 %if no name is assigned, because it is intergenic or whatever
name1='N/A';
end

end