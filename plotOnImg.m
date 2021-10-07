function [] = plotOnImg(Params, ImRef, calPtX, calPtY, Disp, strain, varible)

subplot(1,2,2);
imagesc(repmat(uint8(ImRef),1,1,3));
[x, y] = ndgrid(calPtX,calPtY);
x      = x + reshape(Disp(:,1),Params.Lx,Params.Ly);
y      = y + reshape(Disp(:,2),Params.Lx,Params.Ly);

cLevel = 32;
h1 = subplot(1,2,2);
hold on,
switch varible
    case 'u'
        [C,h1] = contourf(y,x,reshape(Disp(:,1),Params.Lx,Params.Ly),cLevel);
        colorRange = [min(Disp(:,1)),max(Disp(:,1))];
    case 'v'
        [C,h1] = contourf(y,x,reshape(Disp(:,2),Params.Lx,Params.Ly),cLevel);
        colorRange = [min(Disp(:,2)),max(Disp(:,2))];
    case 'exx'
        [C,h1] = contourf(y,x,reshape(strain(:,1),Params.Lx,Params.Ly),cLevel);
        colorRange = [min(strain(:,1)),max(strain(:,1))];
    case 'eyy'
        [C,h1] = contourf(y,x,reshape(strain(:,2),Params.Lx,Params.Ly),cLevel);
    colorRange = [min(strain(:,2)),max(strain(:,2))];
end

colormap('jet')
caxis(colorRange);
colorbar('eastoutside');
set(h1,'LineColor','none');
axis('equal');axis('tight'); 