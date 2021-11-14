clc,clear all
R=100;
L=550e-3;
E=100;
w=2*pi*50;
Z=sqrt(R^2+(w*L)^2);
to=L/R;
alfa=pi/3;
fi=atan(2*pi*50*to);
V=220;
Vm=V*sqrt(2);
syms C;
A=solve(((Vm/Z)*sin(pi/6+alfa-fi)+C*exp(-(pi/6+alfa)/(w*to))-E/R)==((Vm/Z)*sin(5*pi/6+alfa-fi)+C*exp(-(5*pi/6+alfa)/(w*to))-E/R));
Cs=double(A)
sonuc=((Vm/Z)*sin(alfa+pi/6-fi)+Cs*exp(-(alfa+pi/6)/(w*to))-E/R);
if sonuc<0 
        disp('Akım kesintili');
        clear A C Cs;
        syms C
        A=vpasolve(((Vm/Z)*sin(alfa+pi/6-fi)+C*exp(-(alfa+pi/6)/(w*to))-E/R)==0);
        disp('Yeni C değeri bulundu');
        Cs=double(A)
end
if sonuc>0
        
        disp('Akım kesintisiz!')
        
end

teta=fsolve((@(teta)((Vm/Z)*(sin(teta-fi)))+(Cs*exp(-teta/(w*to))-(E/R))),pi)
wt=linspace(alfa+pi/6,teta);
It1=(Vm/Z)*sin(wt-fi)+Cs*exp(-wt/(w*to))-(E/R);
Vdenk= sin(wt)*Vm;
Vyuk_ort=(3/(2*pi))*(trapz(wt,Vdenk)+(E*(alfa+(5*pi/6)-teta)))
IT1_ort=(Vyuk_ort-E)/(3*R)
Iyuk_ort=Vyuk_ort/(3*R)
IT1_etk=sqrt( (1/(2*pi) * trapz(wt,It1.^2) )  )
Iseb_etkin=sqrt(2*IT1_etk.^2)
Pyuk_ort=(1/pi)*trapz(wt,Vdenk.*It1)
