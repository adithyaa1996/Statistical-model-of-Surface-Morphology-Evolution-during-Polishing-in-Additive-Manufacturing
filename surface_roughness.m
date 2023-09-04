function [Heights_for_Sa_calc, mean_model, stdev_model, Sa_model, Hurst_model, KS_test, p_value] = surface_roughness(surface_radii, initial_asperity_radii, h, N, data_stage, packing_density)

Heights_for_Sa_calc = [];

for i = 1:1:N
    Heights_for_Sa_calc = [Heights_for_Sa_calc; -i*ones(round((initial_asperity_radii(i)^2)/20)+2,1)];
end

k = round((1 - packing_density)*sum((initial_asperity_radii.^2))/40);
misc_heights = (quantile(h,0)/2) + ((quantile(h,0.2) - quantile(h,0)).*rand(k,1));
Heights_for_Sa_calc = [Heights_for_Sa_calc; misc_heights];

for i = 1:1:N
    index = find(Heights_for_Sa_calc == -i);
    m = length(index);
    asperity_surface_ratio = (surface_radii(i)^2 / (initial_asperity_radii(i)^2));
    
    if (asperity_surface_ratio < 0.8)
        q1 = round(m*asperity_surface_ratio)+1;
        q2 = round((m - q1)/2);
        Heights_for_Sa_calc(index(1):index(q1)) = h(i);
        if (q2 == 1)
            Heights_for_Sa_calc(index(q1+1): index(end)) = h(i) - ((2/5)*sqrt(initial_asperity_radii(i)^2 - surface_radii(i)^2));
        end
        if (q2 > 1)
            Heights_for_Sa_calc(index(q1+1): index(q2)) = h(i) - ((2/5)*sqrt(initial_asperity_radii(i)^2 - surface_radii(i)^2));
            Heights_for_Sa_calc(index(q2+1): index(end)) = h(i) - (sqrt(initial_asperity_radii(i)^2 - surface_radii(i)^2));
        end     
    end
    
    if (asperity_surface_ratio >= 0.8)
       Heights_for_Sa_calc(index(1):index(end)) = h(i);
    end
end

void_spaces = find(surface_radii > initial_asperity_radii);
active_nodes = length(void_spaces);
void_space_ratio = active_nodes / N;
Heights_for_Sa_calc(end - round(k*void_space_ratio)+1:end) = max(Heights_for_Sa_calc);


mean_model = mean(Heights_for_Sa_calc);
stdev_model = std(Heights_for_Sa_calc);
Hurst_model = Gen_hurst(Heights_for_Sa_calc);

M = length(Heights_for_Sa_calc);
Sa_model = (1/M)*sum(abs(Heights_for_Sa_calc - mean_model),1);

[KS_test, p_value] = KSDivtest(data_stage, max(h) - Heights_for_Sa_calc);