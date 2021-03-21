function w=cholincprec(v)
global P R

w=(v'/R);
w=R\(w');
