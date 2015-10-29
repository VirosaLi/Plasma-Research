function G = Gmatrix( axx,axy,axz,ayx,ayy,ayz,azx,azy,azz )
G(1,:) = [0 0 0 1 0 0];
G(2,:) = [0 0 0 0 1 0];
G(3,:) = [0 0 0 0 0 1];
G(4,:) = [axx axy axz 0 1 0];
G(5,:) = [ayx ayy ayz -1 0 0];
G(6,:) = [azx azy azz 0 0 0];
end

