%% Plotting Workzones
OBJ = ["Wall" "Column"];
Space_BIM=readmatrix('space_BIM.txt');
Err_BIM=readmatrix('err_BIM.txt');
Score_BIM=zeros(length(Err_BIM),1);
az=unique(Space_BIM(:,1));
for i=1:length(Err_BIM)
    x = Space_BIM(Space_BIM(:,1)==az(i),2); % X-coordinates
    y = Space_BIM(Space_BIM(:,1)==az(i),3); % Y-coordinates
    Score_BIM(i,1)=Err_BIM(Err_BIM(:,1)==az(i),2)^2/polyarea(x,y)^1.5*3;
end
x_min = min(Score_BIM);
x_max = max(Score_BIM);
a=8; b=25;
% Avoid division by zero if all values in x are the same
if x_max == x_min
    Score_BIM_C = a * ones(size(Score_BIM)); % If all values are the same, set all to 'a'
else
    Score_BIM_C = a + ((Score_BIM - x_min) / (x_max - x_min)) * (b - a);
end
Score_BIM_C=ceil(Score_BIM_C);

%% Colour map
x_min = min(Score_BIM_C);
x_max = max(Score_BIM_C);
x_normalized = (Score_BIM_C - x_min) / (x_max - x_min);

% Generate a colormap (Jet colormap)
numColors = 256; % Number of colors in the colormap
cmap = jet(numColors);

% Map normalized x values to Jet colormap indices
color_indices = round(x_normalized * (numColors - 1)) + 1;
colors = cmap(color_indices, :); % Extract corresponding RGB values

maxNum=[8,20];
L=sum(maxNum);
stz=Reading_STL_Elements(maxNum,OBJ);

figure; hold
for i=1:L
    [Att_BIM,Bounds_BIM,Area_BIM,BIM_mean]=readSTLPolyGraph(stz(i));
    mx=min(BIM_mean(:,1:2));
    Mx=max(BIM_mean(:,1:2));
    Bx=[mx Mx-mx];
    if isempty(intersect(i,az))
        rr=rectangle('Position',Bx,'Curvature',0.1,'FaceColor',[0 1 0],'EdgeColor','k');
        rr.FaceColor = [0, 1, 0, 0.75];
    else
        x = Space_BIM(Space_BIM(:,1)==i,2); % X-coordinates
        y = Space_BIM(Space_BIM(:,1)==i,3); % Y-coordinates
        p = patch(x, y, colors(az==i,:), 'FaceAlpha', 0.25); % 30% transparent red
        
        % Apply hatching (lines, cross, etc.)
        h = hatchfill2(p, 'single', 'HatchDensity', 20, 'HatchColor', 'k'); %% must be downloaded seperately
        rr=rectangle('Position',Bx,'Curvature',0.1,'FaceColor',[1 0 0],'EdgeColor','k');
        rr.FaceColor = [1, 0, 0, 0.75];
    end
    axis equal;
end
colormap(jet);
colorbar;
clim([x_min x_max]); % Set colorbar scale to original x range

xlabel('x(m)'); ylabel('y(m)'); ylim([-5 15]); xlim([-1 35])
ax = gca; % Get current axes
ax.FontSize = 24; % Set font size
ax.FontName = 'Arial'; % Set font name
