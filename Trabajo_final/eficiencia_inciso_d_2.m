%Con esta funcion se resuelve para el valor de tauy que resulta en
%lss=0.45, manteniendo k/l constante

function F=eficienciamdol1(tauy1,lss2,kss2,A,alpha,delta,gama)

F=((1-alpha)*(1-tauy1)*((A*(kss2^alpha)*(lss2^(1-alpha)))/lss2))-(gama*((A*(kss2^alpha)*(lss2^(1-alpha)))-(kss2-(1-delta)*kss2))/((1-gama)*(1-lss2)));

end
