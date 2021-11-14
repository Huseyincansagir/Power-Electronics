clc,clear all
R=10;
L=10e-3;
w=2*pi*50;
Z=sqrt(R^2+(w*L)^2);
to=L/R;
alfa=pi/2;
fi=atan(2*pi*50*to);
V=220;
Vm=V*sqrt(2)*sqrt(3);
syms C;
A=vpasolve((Vm/Z)*sin(pi/3+alfa-fi)+C*exp(-(pi/3+alfa)/(w*to))==(Vm/Z)*sin(2*pi/3+alfa-fi)+C*exp(-(2*pi/3+alfa)/(w*to)));
Cs=double(A)
sonuc=(Vm/Z)*sin(alfa-fi)+Cs*exp(-alfa/(w*to));
if sonuc<0 
        disp('Akım kesintili');
        clear A C Cs;
        syms C
        A=vpasolve((Vm/Z)*sin(pi/3+alfa-fi)+C*exp(-(pi/3+alfa)/(w*to))==0);
        disp('Yeni C değeri bulundu');
        Cs=double(A)
end
if sonuc>0
        
        disp('Akım kesintisiz!')
        
end

teta=fsolve(@(teta)((Vm/Z)*(sin(teta-fi)))+(Cs*exp(-teta/(w*to))),pi)
wt=linspace(pi/3+alfa,teta);

It1=(Vm/Z)*sin(wt-fi)+Cs*exp(-wt/(w*to));
Vdenk=sin(wt)*Vm;
Vyuk_ort=(3/pi)*trapz(wt,Vdenk)
IT1_ort=Vyuk_ort/(3*R)
Iyuk_ort=Vyuk_ort/2
IT1_etk=sqrt( (1/(pi) * trapz(wt,It1.^2) )  )
Iseb_etkin=sqrt(2*IT1_etk.^2)
Pyuk_ort=(3/pi)*trapz(wt,Vdenk.*It1)
