function [mapa]=mapaRWB()
a = linspace(0,1,128);
g = exp(-((0:255)-128).^2/(2*46^2));

b = [a,ones(1,128)];
r = [ones(1,128),flip(a)];

mapa = [r',g',b']; mapa = flipud(mapa);

