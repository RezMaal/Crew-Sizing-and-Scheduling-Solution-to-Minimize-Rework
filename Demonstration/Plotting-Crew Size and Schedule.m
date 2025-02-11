%% Plotting Crew Sizing and Scheduling Results
figure;
t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax=nexttile;
imagesc(xx_C); % Display matrix as image
colormap(jet); % Use 'jet' colormap
colorbar; % Show color scale
clim(ax,[mX MX]);
title(ax,'Crew Sizing-No Spatial-No Takt');
text_label = sprintf('Duration = %d, Cost = %.2f', t_C, Co_C);
xlabel(ax,text_label)
ax.FontSize = 30; % Set font size
ax.FontName = 'Arial'; % Set font name

ax=nexttile;
imagesc(xx_Cs); % Display matrix as image
colormap(jet); % Use 'jet' colormap
colorbar; % Show color scale
clim(ax,[mX MX]);
title(ax,'Crew Sizing-with Spatial-No Takt');
text_label = sprintf('Duration = %d, Cost = %.2f', t_Cs, Co_Cs);
xlabel(ax,text_label)
ax.FontSize = 30; % Set font size
ax.FontName = 'Arial'; % Set font name

ax=nexttile;
imagesc(xx_Ct); % Display matrix as image
colormap(jet); % Use 'jet' colormap
colorbar; % Show color scale
clim(ax,[mX MX]);
title(ax,'Crew Sizing-No Spatial-with Takt');
text_label = sprintf('Duration = %d, Cost = %.2f', t_Ct, Co_Ct);
xlabel(ax,text_label)
ax.FontSize = 30; % Set font size
ax.FontName = 'Arial'; % Set font name

ax=nexttile;
imagesc(xx_Cts); % Display matrix as image
colormap(jet); % Use 'jet' colormap
colorbar; % Show color scale
clim(ax,[mX MX]);
title(ax,'Crew Sizing-with Spatial-with Takt');
text_label = sprintf('Duration = %d, Cost = %.2f', t_Cts, Co_Cts);
xlabel(ax,text_label)
ax.FontSize = 30; % Set font size
ax.FontName = 'Arial'; % Set font name