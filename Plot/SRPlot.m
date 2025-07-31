% Create plots.
t = tiledlayout(3,3);

min1=min(p_real_grid_pi_6,[],'all');
max1=max(p_real_grid_pi_6,[],'all');
min2=min(p_real_grid_pi_10,[],'all');
max2=max(p_real_grid_pi_10,[],'all');
min3=min(p_real_grid_pi_20,[],'all');
max3=max(p_real_grid_pi_20,[],'all');
min4=min(p1_grid_pi_6,[],'all');
max4=max(p1_grid_pi_6,[],'all');
min5=min(p1_grid_pi_10,[],'all');
max5=max(p1_grid_pi_10,[],'all');
min6=min(p1_grid_pi_20,[],'all');
max6=max(p1_grid_pi_20,[],'all');
min7=min(p1_grid_pi_6_plane,[],'all');
max7=max(p1_grid_pi_6_plane,[],'all');
min8=min(p1_grid_pi_10_plane,[],'all');
max8=max(p1_grid_pi_10_plane,[],'all');
min9=min(p1_grid_pi_20_plane,[],'all');
max9=max(p1_grid_pi_20_plane,[],'all');

min0=min([min1,min2,min3,min4,min5,min6,min7,min8,min9])*0;
max0=max([max1,max2,max3,max4,max5,max6,max7,max8,max9]);

ax1 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p_real_grid_pi_6(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('A slit with a width of','FontSize',25);

ax2 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p_real_grid_pi_10(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('A slit with a width of','FontSize',25);

ax3 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p_real_grid_pi_20(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('A slit with a width of','FontSize',25);

ax4 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p1_grid_pi_6(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
%title('A slit with a width of','FontSize',25);

ax5 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p1_grid_pi_10(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
%title('A slit with a width of','FontSize',25);

ax6 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p1_grid_pi_20(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
%title('A slit with a width of','FontSize',25);

ax7 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p1_grid_pi_6_plane(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
%title('A slit with a width of','FontSize',25);

ax8 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p1_grid_pi_10_plane(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
%title('A slit with a width of','FontSize',25);

ax9 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),p1_grid_pi_20_plane(:,:,3),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
%title('A slit with a width of','FontSize',25);





cb = colorbar;
cb.Layout.Tile = 'south'; 

% Link the axes
t.Padding = 'compact';
t.TileSpacing = 'compact';
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],'xyz');

