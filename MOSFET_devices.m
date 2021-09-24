%%surface voltage%%
qs=-0.5:0.001:1.5;

t=300;
k=1.38*10^-23;
q=1.6*10^-19;
eo=8.854187817*10^-14;
h=6.626*10^-34;
mnmp=1.396*10^-68;

es=11.7;
ni300=1.45*10^10;
Na=6.5*10^13;
eox=3.9;
tox=4.1*10^-7;%%in cm
tox2=25*10^-7;%%in cm
Vgmax=2.5;
Qoi=0;
qms=0;
Mn=600;
w=10*10^-4;%%in cm
L=1*10^-4;%%in cm
vg1=0.525;
vg2=0.5;
vg3=0.475;%%don't use vg less than vth %%you can get vth after simulation the choose vg

ni=2*(2*pi*k*t/h^2)^(3/2) * (mnmp)^(3/4) * exp(-1.12*q/(2*k*t));
qf=(t*k/q)*log(Na/ni);
Cox=eox*eo/tox;
Cox2=eox*eo/tox2;
Vfb=-Qoi*q/Cox+qms;
Vfb2=-Qoi*q/Cox2+qms;
vth=Vfb+sqrt(2*q*es*eo*Na*2*qf)/Cox+2*qf;
vth2=Vfb2+sqrt(2*q*es*eo*Na*2*qf)/Cox2+2*qf;
lamda=(2*10^-6)/L;
%%vds=vg1-vth;

Qs=-(sqrt(2.*es.*eo.*k.*t.*Na).*sqrt((exp(-q.*qs./(k.*t))+q.*qs./(k.*t)-1)+((ni.^2./Na.^2).*(exp(q.*qs./(k.*t))-q.*qs./(k.*t)-1))));
subplot(2,2,1)
semilogy(qs,-Qs);
title('log scale charge - surface potential');
hold on
xline(2*qf,'-.b');
%%xline(vth,'-.r');%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%if U want the graph of -ve voltage remove the "%%" from the else statment
subplot(2,2,2)
m=0;
n=0;
i=1;
for qtemp=qs
    if (qtemp>=0)
        Vg(i)=-Qs(i).*tox./(eox.*eo) + qtemp + Vfb;
    else
        Vg(i)=Qs(i).*tox./(eox.*eo) + qtemp + Vfb;
    end
    if (Vg(i)>(vg1-0.1)&&Vg(i)<(vg1+0.1))
        qs1=qtemp;
    end
    i=i+1;
end
hold on
plot(Vg,qs);
xlim([0 Vgmax])
title('surface potential - gate voltage     with different oxide thickness')
hold on
i=1;
for qtemp=qs
    if (qtemp>=0)
        Vg(i)=-Qs(i).*tox2./(eox.*eo) + qtemp + Vfb2;
    else
        Vg(i)=Qs(i).*tox2./(eox.*eo) + qtemp + Vfb2;
    end
    i=i+1;
end
hold on
plot(Vg,qs);
xlim([0 Vgmax])
ylim([0 inf])
if(Vg(1)==0)
    ylim([0 inf])
end
hold on
yline(2*qf,'-.b');
xline(vth,'-.r');
xline(vth2,'-.r');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)

%%x=sqrt(2*q*es*eo*Na)/Cox;
%%vdmax=(-x/2+sqrt(x^2+4*(vg1-Vfb))/2)^2-2*qf;
vov=vg1-vth;
%%b=vov*(vov*lamda-2);
b=vov*(-2);
c=vov^2;
vdmax=(-b-sqrt(b^2-4*c))/2;
tm=vg1-2*qf-Vfb;
b=-(2*tm+2*q*Na*es*eo/(Cox)^2);
c=tm^2-2*q*Na*es*eo*2*qf/(Cox)^2;
vdmax2=(-b-sqrt(b^2-4*c))/2;
vd=0:0.001:vdmax;
vdd=0:0.001:vdmax2;
Id=Mn.*Cox.*(w/L).*((vg1-Vfb-2*qf-vdd/2).*vdd-(2/3).*(sqrt(2*q*es*eo*Na)/Cox).*((vdd+2*qf).^(3/2)-(2*qf).^(3/2))).* (1+lamda.*vdd);
plot(vdd,Id)
hold on
Id2=Mn.*Cox.*(w/L).*((vg1-vth).*vd-0.5.*vd.^2).* (1+lamda.*vd);
plot(vd,Id2)
title('drain current - drain voltage     with different gate voltage')
%%lamda=Id(end)/(vdmax-Va);
vd2=vdmax:0.001:0.8;
%%lamda=((2.*eo.*es./(q.*Na)).^0.5 ./ (2.*L.*(vd2-vg1+vth+0.7).^2));
Id3=Mn.*Cox.*(w/(2*L)).*(vg1-vth).^2 .* (1+lamda.*vd2);
hold on
plot(vd2,Id3)
vdd2=vdmax2:0.001:0.8;
%%lamda=((2.*eo.*es./(q.*Na)).^0.5 ./ (2.*L.*(vdd2-vg1+vth+0.7).^2));
Id4=Mn.*Cox.*(w/L).*((vg1-Vfb-2*qf-vdmax2/2).*vdmax2-(2/3).*(sqrt(2*q*es*eo*Na)/Cox).*((vdmax2+2*qf).^(3/2)-(2*qf).^(3/2))).* (1+lamda.*vdd2);
hold on
plot(vdd2,Id4)
vds = vdmax;

hold on
%%x=sqrt(2*q*es*eo*Na)/Cox;
%%vdmax=(-x/2+sqrt(x^2+4*(vg2-Vfb))/2)^2-2*qf;
vov=vg2-vth;
b=vov*(-2);
c=vov^2;
vdmax=(-b-sqrt(b^2-4*c))/2;
tm=vg2-2*qf-Vfb;
b=-(2*tm+2*q*Na*es*eo/(Cox)^2);
c=tm^2-2*q*Na*es*eo*2*qf/(Cox)^2;
vdmax2=(-b-sqrt(b^2-4*c))/2;
vd=0:0.001:vdmax;
vdd=0:0.001:vdmax2;
Id=Mn.*Cox.*(w/L).*((vg2-Vfb-2*qf-vdd/2).*vdd-(2/3).*(sqrt(2*q*es*eo*Na)/Cox).*((vdd+2*qf).^(3/2)-(2*qf).^(3/2))).* (1+lamda.*vdd);
plot(vdd,Id)
hold on
Id2=Mn.*Cox.*(w/L).*((vg2-vth).*vd-0.5.*vd.^2).* (1+lamda.*vd);
plot(vd,Id2)
%%lamda=Id(end)/(vdmax-Va);
vd2=vdmax:0.001:0.8;
%%lamda=((2.*eo.*es./(q.*Na)).^0.5 ./ (2.*L.*(vd2-vg2+vth+0.7).^2));
Id3=Mn.*Cox.*(w/(2*L)).*((vg2-vth).^2*(1+lamda.*vd2));
hold on
plot(vd2,Id3)
vdd2=vdmax2:0.001:0.8;
%%lamda=((2.*eo.*es./(q.*Na)).^0.5 ./ (2.*L.*(vdd2-vg2+vth+0.7).^2));
Id4=Mn.*Cox.*(w/L).*((vg2-Vfb-2*qf-vdmax2/2).*vdmax2-(2/3).*(sqrt(2*q*es*eo*Na)/Cox).*((vdmax2+2*qf).^(3/2)-(2*qf).^(3/2))).* (1+lamda.*vdd2);
hold on
plot(vdd2,Id4)

hold on
%%x=sqrt(2*q*es*eo*Na)/Cox;
%%vdmax=(-x/2+sqrt(x^2+4*(vg3-Vfb))/2)^2-2*qf;
vov=vg3-vth;
b=vov*(-2);
c=vov^2;
vdmax=(-b-sqrt(b^2-4*c))/2;
tm=vg3-2*qf-Vfb;
b=-(2*tm+2*q*Na*es*eo/(Cox)^2);
c=tm^2-2*q*Na*es*eo*2*qf/(Cox)^2;
vdmax2=(-b-sqrt(b^2-4*c))/2;
vd=0:0.001:vdmax;
vdd=0:0.001:vdmax2;
Id=Mn.*Cox.*(w/L).*((vg3-Vfb-2*qf-vdd/2).*vdd-(2/3).*(sqrt(2*q*es*eo*Na)/Cox).*((vdd+2*qf).^(3/2)-(2*qf).^(3/2))).* (1+lamda.*vdd);
plot(vdd,Id)
hold on
Id2=Mn.*Cox.*(w/L).*((vg3-vth).*vd-0.5.*vd.^2).* (1+lamda.*vd);
plot(vd,Id2)
%%lamda=Id(end)/(vdmax-Va);
vd2=vdmax:0.001:0.8;
%%lamda=((2.*eo.*es./(q.*Na)).^0.5 ./ (2.*L.*(vd2-vg3+vth+0.7).^2));
Id3=Mn.*Cox.*(w/(2*L)).*((vg3-vth).^2*(1+lamda.*vd2));
hold on
plot(vd2,Id3)
vdd2=vdmax2:0.001:0.8;
%%lamda=((2.*eo.*es./(q.*Na)).^0.5 ./ (2.*L.*(vdd2-vg3+vth+0.7).^2));
Id4=Mn.*Cox.*(w/L).*((vg3-Vfb-2*qf-vdmax2/2).*vdmax2-(2/3).*(sqrt(2*q*es*eo*Na)/Cox).*((vdmax2+2*qf).^(3/2)-(2*qf).^(3/2))).* (1+lamda.*vdd2);
hold on
plot(vdd2,Id4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vgs=vth:0.001:vg1;
Id5=Mn.*Cox.*(w/(2*L)).*((vgs-vth).^2 .* (1+lamda.*vds));
hold on
plot(vgs-vth,Id5)
subplot(2,2,4)
hold on
plot(vgs,Id5)
title('drain current - gate voltage')
