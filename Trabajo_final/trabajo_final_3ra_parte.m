%Este archivo calcula las variables en estado estacionario;
%se usa para calibrar un nivel de gama y de impuestos que den como resultado
%las horas trabajadas en Mexico. Se supone un impuesto a los ingresos
%laborales, regresado a los agentes como transferencia lumpsum.
clear
clc

%Parametros del modelo
A=1;
alpha=0.40;
beta=0.987;
delta=0.012;
gama=0.64;
tauy=0;

%Calculo de variables en estado estacionario
lss=1/(((gama+(1-alpha)*(1-tauy)*(1-gama))*(1-beta*(1-delta))-alpha*beta*delta*gama)/((1-alpha)*(1-gama)*(1-tauy)*(1-beta*(1-delta))));
kss=lss*((A*beta*alpha)/((1-beta*(1-delta))))^(1/(1-alpha));
yss=A*(kss^alpha)*(lss^(1-alpha));
iss=delta*kss;
css=yss-iss;
ratiokl=kss/lss;

%Ahora resolvemos la condicion de eficiencia del mercado laboral para
%obtener el valor de gama que resulta en lss=0.45 (OECD, 2013), manteniendo k/l
%constante:
lss1=0.45;
kss1=lss1*((A*beta*alpha)/((1-beta*(1-delta))))^(1/(1-alpha));
gamaprima=fsolve(@(gamaprima)eficienciamdol(gamaprima,tauy,lss1,kss1,A,alpha,delta),gama);

%Calculo de variables en estado estacionario con nuevo valor de gama
yss1=A*(kss1^alpha)*(lss1^(1-alpha));
iss1=delta*kss1;
css1=yss1-iss1;
ratiokl1=kss1/lss1;

%Ahora resolvemos la condicion de eficiencia del mercado laboral para
%obtener el valor de tauy que resulta en lss=0.45, manteniendo k/l constante:
lss2=lss1;
kss2=kss1;
tauy1=fsolve(@(tauy1)eficienciamdol1(tauy1,lss2,kss2,A,alpha,delta,gama),tauy);

%Calculo de variables en estado estacionario con nuevo valor de tauy
yss2=A*(kss2^alpha)*(lss2^(1-alpha));
iss2=delta*kss1;
css2=yss2-iss2;
ratiokl2=kss1/lss1;

%El hayazgo es consistente con Prescott 2004, a continuacion una curva que
%describe las ofertas de trabajo a distintas tasas impositivas;
grid=100;
lssg=zeros(grid,1);
taus=linspace(-1.4,0.3,grid);

for i=1:grid
    
    lssg(i)=1/(((gama+(1-alpha)*(1-taus(i))*(1-gama))*(1-beta*(1-delta))-alpha*beta*delta*gama)/((1-alpha)*(1-gama)*(1-taus(i))*(1-beta*(1-delta))));
    
end

plot(taus,lssg);
axis([-1.4 0.3 0.226 0.5])
