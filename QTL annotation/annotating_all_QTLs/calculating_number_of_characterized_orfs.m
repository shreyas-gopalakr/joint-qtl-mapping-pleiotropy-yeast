V = GFFAnnotation('BY4742_sorted_annotations.gff3');
orf_classifications = {};

for j = 1:V.NumEntries
    gene_feature = V.Feature(j);
    if strcmp(gene_feature{1,1}, 'gene') == 1
        attr = V.Attributes(j);
        attr_str = attr{1,1};
        attr_split = split(attr_str, ';');
        
        % search for orf_classification= rather than assuming position
        orf_class = '';
        for k = 1:length(attr_split)
            if contains(attr_split{k,1}, 'orf_classification=')
                orf_temp = split(attr_split{k,1}, 'orf_classification=');
                orf_class = orf_temp{2,1};
                break
            end
        end
        
        if isempty(orf_class)
            orf_class = 'not found';
        end
        
        orf_classifications{end+1} = orf_class;
    end
end

% count occurrences of each unique classification
[unique_classes, ~, ic] = unique(orf_classifications);
for k = 1:length(unique_classes)
    fprintf('%s: %d\n', unique_classes{k}, sum(ic == k));
end