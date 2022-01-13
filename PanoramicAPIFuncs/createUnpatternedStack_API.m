function z0 = createUnpatternedStack_API(stack,coordsX,coordsY)
    z0 = 0;
    %% Repeating layers
    for i = 1:stack.Periodic.N
        for iMat = length(stack.Periodic.material):-1:1 % Reverse order, because Panoramic wants bottom to top
            %% Add layer
            t = stack.Periodic.t(iMat)*stack.Periodic.dspace;
            if t > 0
                coords = [coordsX coordsY z0 z0+t];
                z0 = z0+t;
    %             fprintf(['Layer ' num2str(i) ', Mat ' num2str(iMat) ', z = ' num2str(z0) '\n'])
                coordString = coords2String(coords);
                material = stack.Periodic.material{iMat};
                callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);
            end
        end
    end

    %% Capping layers
    for iMat = length(stack.Cap.material):-1:1 % Reverse order, because Panoramic wants bottom to top
        %% Add layer
        t = stack.Cap.t(iMat);
        if t > 0
            coords = [coordsX coordsY z0 z0+t];
            z0 = z0+t;
            coordString = coords2String(coords);
            material = stack.Cap.material{iMat};
            callJavaSubroutine('panoramictech.v700.OpenAPI.CAPIClient','addBlock','box',material,coordString);
        end
    end
end