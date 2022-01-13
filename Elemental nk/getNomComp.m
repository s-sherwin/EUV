function composition = getNomComp(layer_mat,materials)
    composition = zeros(length(layer_mat),length(materials)); % Nominal composition based on stoichiometry
    for i = 1:size(composition,1)
        str = layer_mat{i};
%         digits = isstrprop(str,'digit');
        uppers = isstrprop(str,'upper');
        uppers = find(uppers);
        for j = 1:length(uppers)
            i0 = uppers(j);
            if j == length(uppers)
                i1 = length(str);
            else
                i1 = uppers(j+1)-1;
            end
            sub_str = str(i0:i1);
            digits = isstrprop(sub_str,'digit');
            mat_ind = strcmp(materials,sub_str(~digits));
            if any(digits)
                mat_num = str2double(sub_str(digits));
            else
                mat_num = 1;
            end
            composition(i,mat_ind) = mat_num;
        end
    end
    composition = composition./sum(composition,2);
    composition = composition';
end