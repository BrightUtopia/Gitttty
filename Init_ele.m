%%%生成https://doi.org/10.1016/j.cma.2020.113307 相似的初始构型

function xPhys=Init_ele(nelx,volfrac)

xx=0.5:1:nelx;
yy=xx;
[ye,xe]=meshgrid(xx,yy);
parfor i=1:nelx
    for j=1:nelx
        xP(i,j)=1/2*(sin(pi*xe(i,j)/nelx)+sin(pi*ye(i,j)/nelx));
    end
end

c=sum(sum(xP))/(volfrac*nelx*nelx);
xPhys=xP./c;
% figure
% colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; drawnow;
% xlim([0.5,nelx+0.5]);xticks([]);yticks([]);

end