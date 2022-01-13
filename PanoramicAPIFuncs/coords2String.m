function coordString = coords2String(coords)
    coordString = '';
    for i = 1:length(coords)
        coordString = [coordString num2str(coords(i))];
        if i < length(coords)
            coordString = [coordString ','];
        end
    end
end