function [phi,theta] = cart2torus(x,R,r)

     [phi,~] = cart2pol((sqrt(x(:,1).^2+x(:,2).^2)-R)/r,x(:,3)/r);
     
     [theta,~] = cart2pol(x(:,1)./(R+r*cos(phi)),x(:,2)./(R+r*cos(phi)));

end