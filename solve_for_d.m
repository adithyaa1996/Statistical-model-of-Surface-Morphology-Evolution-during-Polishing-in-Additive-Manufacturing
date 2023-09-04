%% solve for distance
function d_crit = solve_for_d(Force,radii,heights,params,vectorlength)

E = params{1};
dd = params{2};

%  d_init = 0;
%  d = d_init;

gap = zeros(length(vectorlength),1); 
k = 1; 
k_crit = 0;

for d = 0:dd:max(heights)
       active_nodes = find(heights>d);
       gap(k) = Force - ((2/3)*E*sum(sqrt(radii(active_nodes)).*(abs(heights(active_nodes) - d.*ones(length(active_nodes),1))).^(3/2)));

       if (gap(k)>0 && k_crit == 0)
            k_crit = k;
       end
            k = k+1;
end
% end
d_crit = (k_crit-1) * dd-dd/2;
active_nodes = find(heights>d_crit);
% Force - 2/3*sum(sqrt(radii(active_nodes)).*(heights(active_nodes)-d_crit*ones(length(active_nodes),1)).^(3/2))


