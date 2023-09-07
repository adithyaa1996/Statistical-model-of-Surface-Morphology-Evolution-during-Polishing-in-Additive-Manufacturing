%% Initial Surface Generation

% Consider a 0.25mm x 0.25mm surface for the purpose of modelling.
% Surface packing density is assumed to lie between 0.7 and 0.9 based on
% the specified minimun and maximum radii value of asperities.

% Placing asperities on a 2D space
grid_length = 250;  % Units in micrometers
[Locations,radii_m0] = random_circle_packing_rectangle([grid_length grid_length],3.14,25,true);
x_pos = Locations(:,1);
y_pos = Locations(:,2);

figure;
DT = delaunayTriangulation(x_pos, y_pos);
triplot(DT);
title('Surface Asperities - Network Formation');
xlabel('X (10^{-6})m')
ylabel('Y (10^{-6})m')
pub_fig;
