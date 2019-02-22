function z = positive_ellipsoid(a,b,c,x0,y0,x,y)

Iin = (x-x0).^2./a^2 + (y-y0).^2./b^2 <= 1;

z = x-x; % lazy init

z       = sqrt((((x-x0).^2./a^2+(y-y0).^2./b^2) - 1) * -c^2);
z(~Iin) = 0;

end
