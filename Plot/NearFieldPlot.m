% Create plots.
t = tiledlayout(2,3);

min1=min(real(E1k3_near_grid(:,:,3)),[],'all');
max1=max(real(E1k3_near_grid(:,:,3)),[],'all');
min2=min(real(E2k3_near_grid(:,:,3)),[],'all');
max2=max(real(E2k3_near_grid(:,:,3)),[],'all');
min3=min(real(E3k3_near_grid(:,:,3)),[],'all');
max3=max(real(E3k3_near_grid(:,:,3)),[],'all');
min4=min(imag(E1k3_near_grid(:,:,3)),[],'all');
max4=max(imag(E1k3_near_grid(:,:,3)),[],'all');
min5=min(imag(E2k3_near_grid(:,:,3)),[],'all');
max5=max(imag(E2k3_near_grid(:,:,3)),[],'all');
min6=min(imag(E3k3_near_grid(:,:,3)),[],'all');
max6=max(imag(E3k3_near_grid(:,:,3)),[],'all');
min0=min([min1,min2,min3,min4,min5,min6])+0.1;
max0=max([max1,max2,max3,max4,max5,max6])-0.1;

ax1 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),real(E1k3_near_grid(:,:,3)),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('Real part of E1','FontSize',25);

ax2 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),real(E2k3_near_grid(:,:,3)),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('Real part of E2','FontSize',25);

ax3 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),real(E3k3_near_grid(:,:,3)),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('Real part of E3','FontSize',25);

ax4 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),imag(E1k3_near_grid(:,:,3)),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('Imaginary part of E1','FontSize',25);

ax5 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),imag(E2k3_near_grid(:,:,3)),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('Imaginary part of E2','FontSize',25);

ax6 = nexttile;
mesh(X_near_grid(:,:,1),Y_near_grid(:,:,1),imag(E3k3_near_grid(:,:,3)),'FaceColor','flat');
set(gca,'XColor', 'none','YColor','none');
clim([min0,max0]);
title('Imaginary part of E3','FontSize',25);

cb = colorbar;
cb.Layout.Tile = 'south'; 

% Link the axes
t.Padding = 'compact';
t.TileSpacing = 'compact';
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'xyz');

