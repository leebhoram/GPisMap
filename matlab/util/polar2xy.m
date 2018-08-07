function [x,y] = polar2xy(theta,r)

theta = reshape(theta,1,[]);
r = reshape(r,1,[]);
Ct = cos(theta);
St = sin(theta);
x = r.*Ct;
y = r.*St;

if nargout == 1
    x = [x;y];
end

end