function [xp,yp,zp] = mapc2m(xc,yc)

global map isflat;

map = 'cart';
% map = 'torus';
% map = 'latlong';
map = 'annulus';

R = 1;
r = 0.4;

switch map
    case 'cart'
        isflat = true;        
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
                [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_cart(xc1,yc1);
        xp = 5*(xp+1);
        yp = 3*(yp + 1);
        
        brick_data = load('brick.dat');
        mi = brick_data(1,1);
        mj = brick_data(1,2);
        xv = brick_data(2:end,1);
        yv = brick_data(2:end,2);

        blockno = getblocknumber();
        xs = (xv(blockno+1)+1-mi/2)*0.2;
        ys = (yv(blockno+1)+1-mj/2)*0.2;
        xp = xp + xs;
        yp = yp + ys;

    case 'torus'
        isflat = false;
        % Find center to expand radius around.                        
        [xc1,yc1,~] = mapc2m_brick(0.5,0.5);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
        th = atan2(yp,xp);
        phi = atan2(zp,xp-1);        
        v = [xp-cos(th),yp-sin(th),zp];
        s = 0.2;
        
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_torus(xc1,yc1);
        
        % 
        xp = xp + s*(cos(th) + v(1));
        yp = yp + s*(sin(th) + v(2));
        zp = zp + s*v(3);
        
    case 'latlong'
        isflat = false;
        lat = [-50 50];
        a = 0.5+lat(1)/180;
        b = 0.5+lat(2)/180;
        
        [xc1,yc1,~] = mapc2m_brick(0.5,0.5);
        [xp,yp,zp] = mapc2m_latlong(xc1,yc1);
        v = [xp,yp,zp];
        
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        yc1 = a + (b-a)*yc1;        
        [xp,yp,zp] = mapc2m_latlong(xc1,yc1);
        s = 0.2;
        xp = xp + s*v(1);
        yp = yp + s*v(2);
        zp = zp + s*v(3);
    case 'annulus'
        isflat = true;
        [xc1,yc1,~] = mapc2m_brick(0.5,0.5);
        [xp,yp,zp] = mapc2m_annulus(xc1,yc1);
        v = [xp,yp];
        
        [xc1,yc1,~] = mapc2m_brick(xc,yc);
        [xp,yp,zp] = mapc2m_annulus(xc1,yc1);
        s = 0.2;
        xp = xp + s*v(1);
        yp = yp + s*v(2);
end

end