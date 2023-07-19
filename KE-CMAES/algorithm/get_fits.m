function [data,fits,nets,layers]=get_fits(data,nets,layers,pop_I,pro)
    if pro.env+1 <6
        fits_I = pro.GetFits(pop_I);
        fits = fits_I;
        data = [data,fits];
    else
        if pro.env +1 ==6
            % 初始预测
            fits_I = pro.GetFits(pop_I);
            fits = fits_I;
            data = [data,fits];
            options = trainingOptions('adam', 'MaxEpochs', 20, 'MiniBatchSize', 64);
            nets = trainNetwork(data(:,1:5)', data(:,6)', layers, options);
            layers = nets.Layers;
    
        else
            fits_I = pro.GetFits(pop_I(1:2000,:));
%               fits = fits_I;
            [results,nets] = net_train(nets,data,fits_I,pop_I,pro);
%             results = fits_I';
            data = [data,results'];
            fits = results';
        end
    
    end
end