%Con esta funcion se resuelve para el valor de gama que resulta en lss=0.45
%(OECD, 2013), manteniendo k/l constante

function F=eficienciamdol(gamaprima,tauy,lss1,kss1,A,alpha,delta)

F=((1-alpha)*(1-tauy)*((A*(kss1^alpha)*(lss1^(1-alpha)))/lss1))-(gamaprima*((A*(kss1^alpha)*(lss1^(1-alpha)))-(kss1-(1-delta)*kss1))/((1-gamaprima)*(1-lss1)));

end
