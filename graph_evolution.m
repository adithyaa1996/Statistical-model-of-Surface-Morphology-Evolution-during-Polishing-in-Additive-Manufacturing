function [Adj, Fiedler, eigenvector, bins_connectedcomp] = graph_evolution(X, Y, packing_density, Force, r, h, heights_m0, params, vectorlength, stagenumber)

x_pos_new = X;
y_pos_new = Y;
N = length(X);

d_critical = solve_for_d(Force, r, h, params, vectorlength);
[x_active,y_active] = find(heights_m0 > d_critical);
m = length(x_active);
x_index = zeros(m,1);
y_index = zeros(m,1);
for i = 1:length(x_active)
    x_index(i) = x_active(i);
    x_active(i) = x_pos_new(x_index(i));
    y_index(i) = find(x_pos_new == x_active(i));
    y_active(i) = y_pos_new(y_index(i));
end

figure
DT = delaunayTriangulation(x_pos_new,y_pos_new);
triplot(DT,'--');
pub_fig;
set(gca, 'color', [0.8 0.8 0.8])
str = sprintf("Network model: Stage - " + stagenumber);
title(str)
xlabel('x')
ylabel('y')
hold on
% scatter(x_pos_new, y_pos_new, 'green' , 'filled')
% hold on
scatter(x_active, y_active, 20, 'red', 'filled')
hold on

for k = 1:1:length(x_index)
    i = x_index(k);
        for j = 2:1:N
            TF = isConnected(DT,i,j);
            cc_dist = 0;
            sum_of_radii = 0;
            rho = 0;
            
            if (TF == 1)
                cc_dist = sqrt((X(i) - X(j))^2 + (Y(i) - Y(j))^2);
                sum_of_radii = r(i) + r(j);
                rho = max(cc_dist - sum_of_radii, 0);
                edge_prob(i,j) = exp(-pi()*packing_density*rho^2);
                cbar = colorbar;
                set(cbar, 'ylim', [0 1]);
                % Colouring
                if (edge_prob(i,j) < 0.1)
                    plot([X(i),X(j)], [Y(i),Y(j)],'Linewidth', 0.65, 'Color',[0.3, 0.2, 0.8]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.1  && edge_prob(i,j) < 0.2)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 0.65, 'Color', [0.3, 0.3, 0.9]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.2  && edge_prob(i,j) < 0.3)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 0.65, 'Color', [0.2, 0.5, 1]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.3  && edge_prob(i,j) < 0.35)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 0.65, 'Color', [0.2, 0.6, 0.9]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.35  && edge_prob(i,j) < 0.4)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 0.65, 'Color', [0.1, 0.6, 0.9]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.4  && edge_prob(i,j) < 0.45)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 0.65, 'Color', [0.1, 0.7, 0.9]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.45  && edge_prob(i,j) < 0.5)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 0.65, 'Color', [0, 0.7, 0.8]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.5  && edge_prob(i,j) < 0.55)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [0.1, 0.8, 0.7]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.55  && edge_prob(i,j) < 0.6)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [0.2, 0.8, 0.6]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.6  && edge_prob(i,j) < 0.65)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [0.3, 0.8, 0.5]);
                    hold on
                    
                 elseif (edge_prob(i,j) >= 0.65  && edge_prob(i,j) < 0.7)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [0.5, 0.8, 0.3]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.7  && edge_prob(i,j) < 0.75)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [0.7, 0.8, 0.2]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.75  && edge_prob(i,j) < 0.8)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [0.8, 0.7, 0.2]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.8  && edge_prob(i,j) < 0.85)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [1, 0.7, 0.2]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.85  && edge_prob(i,j) < 0.9)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [1, 0.8, 0.2]);
                    hold on
                    
                elseif (edge_prob(i,j) >= 0.9  && edge_prob(i,j) < 0.95)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [1, 0.9, 0.2]);
                    hold on
                
                elseif (edge_prob(i,j) >= 0.95  && edge_prob(i,j) <= 1)
                    plot([X(i),X(j)], [Y(i),Y(j)], 'Linewidth', 1.5, 'Color', [1, 0.9, 0.1]);
                    hold on
                    
                end
            end
        end
end

    % Delaunay triangulation
    tri = delaunay(X,Y);
    % Calculate adjacency matrix
    AdjMat = zeros(N,N);
    for kk = 1:size(tri,1)
        AdjMat(tri(kk,1), tri(kk,2)) = true;
        AdjMat(tri(kk,2), tri(kk,3)) = true;
        AdjMat(tri(kk,3), tri(kk,1)) = true;
    end
    AdjMat = AdjMat | AdjMat';
    
    for i=1:1:size(AdjMat,1)
        for j = 1:1:size(AdjMat,2)
            if (AdjMat(i,j) == 1)
                cc_dist = sqrt((X(i) - X(j))^2 + (Y(i) - Y(j))^2);
                sum_of_radii = r(i) + r(j);
                rho = max(cc_dist - sum_of_radii, 0);
                if (exp(-pi()*packing_density*rho^2) > 0.85)
                    AdjMat(i,j) = 1;
                else
                    AdjMat(i,j) = 0;
                end
            end
        end
    end
    
    G = graph(AdjMat);
    [bins,binsize] = conncomp(G);
    largest_comp = mode(bins);
    idx = find(bins == largest_comp);
    H = subgraph(G, idx);
    
    Adj = adjacency(H);
    D = diag(sum(Adj,2));             % Degree Matrix
    L = D - Adj;                      % Laplacian Matrix             
    L_norm = (D^-0.5)*L*(D^-0.5);      % Normalized L
    E = sort(eig(L_norm));                    % Eigen values of L
    Fiedler = E(2);
    if Fiedler < 0
        Fiedler = 0;
    end
    eigenvector = E;
    
    binsize_new = [];
    
    for i = 1: length(binsize)
        if binsize(i) > 1
            binsize_new(1, end+1) = binsize(i);
            % binsize_new(2, end) = i;
        end
    end
    
    bins_connectedcomp = binsize_new;
    Adj = full(Adj);
                
                