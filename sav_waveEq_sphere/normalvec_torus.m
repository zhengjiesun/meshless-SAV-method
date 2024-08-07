function nr = normalvec_torus(x,R,r)

% compute the normal vetor of torus, with two radius R, r
torus = @(la,th) [(R+r*cos(la))*cos(th) (R+r*cos(la))*sin(th) r*sin(la)];
syms la th;
nr = inline(cross(diff(torus(la,th),la),diff(torus(la,th),th)));
[la, th] = cart2torus([x(:,1),x(:,2),x(:,3)],R,r);
nr = nr(la,th); nr = nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);