function w=ssorprec(v)
global D DomL

om=1.3;
w=D*(DomL\v);
w=(2-om)*(DomL'\w);
