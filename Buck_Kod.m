T=100;
f=1/T;
Vimin=10;
Vimax=20;
Vomin=50;
Vomax=100;
Dmin=1-(Vimax/Vomin);
Dmax=1-(Vimin/Vomax);
R=100;
L=2.4e-3;
P=55;
Iyuk=sqrt(P/R);
dImin=((Vimax-Vomin)*Dmin)/(f*L)
dImax=((Vimax-Vomin)*Dmax)/(f*L)
ILmin1=Iyuk-dImin/2
ILmax1=Iyuk+dImin/2
ILmin2=Iyuk-dImax/2
ILmax2=Iyuk+dImax/2
Lkrit=(1-Dmin)*R/(2*f);
Ckrit=(1-Dmin)/(16*L*f^2);
C=(Vimin*(1-Dmin)*Dmin)/(8*L*0.1*f^2);