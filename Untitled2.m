Newdata=xlsread('1','sheet1');   
Sb=117.8;%密相区通风截面积
Sv=229.15;% m2
h1=4.78;%密相区压力测点之间距离 m
h2=4.8;%稀相区测点至炉膛高度距离 m
Fg_m=(Newdata(:,2)+Newdata(:,8)+Newdata(:,3))/3600;%密相区风量Nm3/s  
Fg_x=(Newdata(:,1))/3600;%稀相区上二次风 Nm3/s
Fg=Newdata(:,11)/3600;%总风量 Nm3/s
Tg1=Newdata(:,23);%空预器出口一次风温度
T1=Newdata(:,9);%密相区床温
T2=Newdata(:,41);%稀相区床温
T3=Newdata(:,69);%分离器温度
% Fg_m1=Fg_m./((273./(273+T1)).*((Newdata(:,43)*1000+101325)/101325));
% Fg1=Fg./((273./(273+T2)).*((150+101325)/101325));
V1=Fg_m.*(273+T1)/(273*Sb);%密相区实际流化速度m/s
V2=Fg.*(273+T2)/(273*Sv);%稀相区实际流化速度m/s
oxy_air=Newdata(:,68)/100;%空预器入口含氧量“%”
dc1=0.002;row_c1=1800;%密相区颗粒直径以及颗粒密度
dc2=0.0002;row_c2=1800;
Vb=1795.5;
Vs=4033;%各相区体积m3a
delta_p1=(Newdata(:,43)-Newdata(:,44))*1000;%密相区床压差，根据测点图确定
delta_p2=Newdata(:,15);%炉膛上部差压
M1=delta_p1*Vb/(9.8*h1);
M2=delta_p2*Vs/(9.8*h2);%各相区床料质量

row_g=1.293;%空气密度在0℃，标准大气压下
row_gm1=row_g*(273./(273+T1)).*((Newdata(:,43)*1000+101325)/101325);%在密相区的温度压力下，空气密度
row_gm2=row_g*(273./(273+T2)).*((150+101325)/101325);
Ks1=3*0.513*(T1+273).*exp(-9160./(273+T1));%燃烧速率Rc1常数
Ks2=3*0.513*(T2+273).*exp(-9160./(T2+273));%燃烧速率Rc2常数
Ks3=0.513*(T3+273).*exp(-9160./(T3+273));
M=(Newdata(:,6)+Newdata(:,7))*1000/3600;%单位时间床料总量kg/s
M_c=Newdata(:,7)*1000/3600;%单位时间内给煤量 kg/s
C_bed=(Newdata(:,7))*1000*0.4103/3600;%燃料量中碳含量 kg/s
%F_hff=(Newdata(:,7))*173.25/3600;%取自空干基混煤挥发分平均值 kg/s
N_fl=Newdata(:,17);%分离效率
Vmf1=0.294*(0.002^0.584./(1.525e-4).^0.056).*((row_c1-row_gm1)./row_gm1).^0.528;%临界流化速度 row_gm1
Vmf2=0.294*(0.0002^0.584./(1.525e-4).^0.056).*((row_c2-row_gm2)./row_gm2).^0.528;
Fg_gas=Newdata(:,21).*Newdata(:,20)/(100*3600);%烟气含氧量 Nm3/s

%挥发份计算
f_hff=-exp(26.41-3.961*log(T1)+1.15*17.325/100)+79*17.325/100+2.289;
V0=M_c*(1-0.01671-0.44408).*f_hff/100;%挥发分总量kg/s
V_CH4=(0.201-0.469*0.17325+0.241*0.17325^2)*V0;
V_H2=(0.157-0.868*0.17325+1.388*0.17325^2)*V0;
V_CO=(0.428-2.653*0.17325+4.845*0.17325^2)*V0;
V_H2O=(0.409-2.389*0.17325+4.554*0.17325^2)*V0;
V_CO2=(0.135-0.9*0.17325+1.906*0.17325^2)*V0;
V_oil=V0-(V_H2+V_CH4+V_CO+V_H2O+V_CO2);%0.325-7.279*0.17325+12.88*0.17325^2;%挥发份中各组分质量份额 CH0.689 O 0.014
Fhff_O2=2*(V_CH4/16)+0.5*(V_H2/2)+0.5*(V_CO/28)+1.165*(V_oil/12.913);%挥发份消耗氧气 kmol/s
Q_hff=V_CH4*50016+V_H2*124238+V_CO*10077+V_oil*37000;%挥发分放出热量KJ

%烟气成分计算
Scon=M_c*0.3722/100;%燃料含S量kg/s
Smoke_s=Newdata(:,21).*Newdata(:,18)/(1000000*3600);%烟气含SO2量 kg/s
desulfu=(Scon-Smoke_s*32/64)./Scon;%脱硫效率 矩阵
Smoke_o2=Fg_gas*32/22.4;%氧气单位kg/s
%Smoke_NOx=Newdata(:,21).*Newdata(:,19)/(1000000*3600);%烟气含NOx量kg/s
Smoke_N=M_c*0.7171*22.4/(100*28)+Newdata(:,11)*0.78/3600;%烟气含N2量Nm3/s
Smoke_W=1.29*Newdata(:,11)*7.23/(1000*3600)+M_c*18*1.38/100+M_c*7.51/100;%烟气含水量kg/s
Smoke_CO2=Newdata(:,21)/3600-Smoke_s*22.4/64-Smoke_o2*22.4/32-Smoke_N-Smoke_W*22.4/18;%烟气CO2含量Nm3/s
row_co2=((101475*273)./(101325*(273+T2)))*1.96;
row_N2=((101475*273)./(101325*(273+T2)))*1.25;
row_so2=((101475*273)./(101325*(273+T2)))*2.86;
row_h2o=((101475*273)./(101325*(273+T2)))*0.8;
row_o2=((101475*273)./(101325*(273+T2)))*1.429;
row_gas=row_co2.*(Smoke_CO2)./(Newdata(:,21)/3600)+row_N2.*(Smoke_N)./(Newdata(:,21)/3600)+row_so2.*((Smoke_s)*22.4)./(64*Newdata(:,21)/3600)+row_h2o.*((Smoke_W)*22.4)./(18*Newdata(:,21)/3600)+row_o2.*((Smoke_o2)*22.4)./(32*Newdata(:,21)/3600);%烟气密度 kg/m3
Row_gas=1.96*(Smoke_CO2)./(Newdata(:,21)/3600)+1.25*(Smoke_N)./(Newdata(:,21)/3600)+2.86*((Smoke_s)*22.4)./(64*Newdata(:,21)/3600)+0.8*((Smoke_W)*22.4)./(18*Newdata(:,21)/3600)+1.429*((Smoke_o2)*22.4)./(32*Newdata(:,21)/3600);%标况下烟气密度
Fg_c=Newdata(:,21).*Row_gas*0.17/(100*3600);%烟气含碳量kg/s

Y=zeros(140,5);
% row_o2=((101475*273)./(101325*(273+T2)))*1.429;%氧气在稀相区的密度
C_o1=0.00006;
C_o2=(0.21*mean(Fg(332:471))/22.4-mean(0.93*C_bed(332:471))/12-mean(oxy_air(332:471).*Newdata(332:471,21)/(3600*22.4))-mean(Fhff_O2(332:471))- C_o1*Vb)/Vs;
%C_o2=(0.21*mean(Fg(577:716))/22.4-mean(0.93*C_bed(577:716))/12-mean(oxy_air(577:716).*Newdata(577:716,21)/(3600*22.4))-mean(Fhff_O2(577:716))- C_o1*Vb)/Vs;
%C_o2=mean(oxy_air/22.4);%稀相区氧气摩尔浓度kmol/m3
%F_out1=C_o2*Vs+mean(oxy_air(577:716).*Newdata(577:716,21)/(3600*22.4))+mean(R_c2(332:812))/12-0.21*mean(Fg_x(332:812))/22.4;%由密相区进入稀相区的氧气总量kmol
% C_o11=(mean(0.23*Fg_m(332:812))/22.4-mean(R_c1(332:812))/12-F_out1-Fhff_O2)./Vb;%密相区迭代浓度 kmol/m3
% while C_o1>0
%     C_o1=C_o1-0.0001
% % while abs(C_o11-C_o1)>0.001
%     
% %   C_o1=(C_o1+0.001*rand(1));
% %     if C_o11>C_o1
% %     C_o1=C_o1-(C_o11-C_o1)/2;
% %     else
% %       C_o1=C_o1+(C_o1-C_o11)/2;  
% %     end
    
for i=332:1:471
   x0=zeros(1,5);
   a1=M(i);
   a2=N_fl(i)*M2(i)*(V2(i)-Vmf2(i));
   a3=72*Ks1(i)*C_o1/(dc1*row_c1);
   a4=M1(i)*(V1(i)-Vmf1(i));
   b1=72*Ks2(i)*C_o2/(dc2*row_c2);
   b2=M2(i)*(V2(i)-Vmf2(i));
   c1=C_bed(i);
   c21=N_fl(i);
   %c2=N_fl(i)*(V2(i)-Vmf2(i));
   c3=V1(i)-Vmf1(i);
   d1=V2(i)-Vmf2(i);
   d2=0.3*Ks3(i);
   d3=Fg_c(i);
  %OPTIONS=optimoptions('fsolve','Algorithm','Levenberg-Marquardt')
   x=fsolve(@(x) myfun2(x,a1,a2,a3,a4,b1,b2,c1,c21,c3,d1,d2,d3),x0);
    Y(i-331,:)=x;
%  a2*x(4)
   a4*x(3)
   b2*x(4)
   b1*x(2)
end
Y
% I=[332:1:471];
% plot(I,Y(:,1))

R_c1=72*C_o1*Ks1(332:471).*Y(:,1)/(dc1*row_c1);%燃烧速率kg/s
R_c2=72*C_o2*Ks2(332:471).*Y(:,2)/(dc2*row_c2);
R_c3=0.3*Ks3(332:471).*Y(:,5);
a=M1.*(V1-Vmf1);
b=M2.*(V2-Vmf2);
Fu1=Y(:,3).*a(332:471);
Fu2=Y(:,4).*b(332:471);
% I=[332:1:471];
% figure(1)
% plot(I,Ks1(332:471),'k')
y=mean(Y);
m_c1=y(1);
m_c2=y(2);
k1=y(3);
k2=y(4);
m_c3=y(5);
R_c1=72*C_o1*Ks1.*m_c1/(dc1*row_c1);%燃烧速率kg/s
R_c2=72*C_o2*Ks2.*m_c2/(dc2*row_c2);
R_c3=0.3*Ks3.*m_c3;
% F_out1=C_o2*Vs+mean(oxy_air(332:471).*Newdata(332:471,21)/(3600*22.4))+mean(R_c2(332:471))/12-0.21*mean(Fg_x(332:471))/22.4;%由密相区进入稀相区的氧气总量kmol
% C_o11=(mean(0.21*Fg_m(332:471))/22.4-mean(R_c1(332:471))/12-F_out1-Fhff_O2)./Vb;%密相区迭代浓度 kmol/m3

Fu2=k2.*M2.*(V2-Vmf2);
Fu1=k1.*M1.*(V1-Vmf1);%密相区进入稀相区的物料量kg/s
% Fcc=30*M_c;
% Fu2=Fcc./N_fl;
% Fu1=Fu2+R_c2;
Fcc=N_fl.*Fu2;%循环物料量
%Uv_o2=0.21*Newdata(:,11)/3600-((R_c1+R_c2)*22.4/12)-(Newdata(:,21).*Newdata(:,20)/3600);%挥发份消耗氧量m3/s
% sum(Newdata(332:812,3)*0.21/3600)
% %sum((Newdata(332:812,1)+Newdata(332:812,2))*0.21/3600)
% sum(R_c1*22.4/12)

%  F_out1=C_o2*Vs+mean(oxy_air(332:812).*Newdata(332:812,21)/(3600*22.4))+mean(R_c2(332:812))/12-0.23*mean(Fg_x(332:812))/22.4;%由密相区进入稀相区的氧气总量kmol
%  C_o11=(mean(0.23*Fg_m(332:812))/22.4-mean(R_c1(332:812))/12-F_out1-Fhff_O2)./Vb;%密相区迭代浓度 kmol/m3
  
%  result(j,:)=[C_o11,C_o1]
%  j=j+1;
% end
% T_gas=Newdata(:,40);%烟气温度℃
% Scon=M_c*0.3722/100;%燃料含S量kg/s
% Smoke_s=Newdata(:,21).*Newdata(:,18)/(1000000*3600);%烟气含SO2量 kg/s
% desulfu=(Scon-Smoke_s*32/64)./Scon;%脱硫效率 矩阵
% Smoke_o2=Fg_gas*32/22.4;%氧气单位kg/s
% %Smoke_NOx=Newdata(:,21).*Newdata(:,19)/(1000000*3600);%烟气含NOx量kg/s
% Smoke_N=M_c*0.7171*22.4/(100*28)+Newdata(:,11)*0.78/3600;%烟气含N2量Nm3/s
% Smoke_W=1.29*Newdata(:,11)*7.23/(1000*3600)+M_c*18*1.38/100+M_c*7.51/100;%烟气含水量kg/s
% Smoke_CO2=Newdata(:,21)/3600-Smoke_s*22.4/64-Smoke_o2*22.4/32-Smoke_N-Smoke_W*22.4/18;%烟气CO2含量Nm3/s
% 
% row_co2=((101475*273)./(101325*(273+T2)))*1.96;
% row_N2=((101475*273)./(101325*(273+T2)))*1.25;
% row_so2=((101475*273)./(101325*(273+T2)))*2.86;
% row_h2o=((101475*273)./(101325*(273+T2)))*0.8;
% row_o2=((101475*273)./(101325*(273+T2)))*1.429;
% row_gas=row_co2.*(Smoke_CO2)./(Newdata(:,21)/3600)+row_N2.*(Smoke_N)./(Newdata(:,21)/3600)+row_so2.*((Smoke_s)*22.4)./(64*Newdata(:,21)/3600)+row_h2o.*((Smoke_W)*22.4)./(18*Newdata(:,21)/3600)+row_o2.*((Smoke_o2)*22.4)./(32*Newdata(:,21)/3600);%烟气密度 kg/m3
% Row_gas=1.96*(Smoke_CO2)./(Newdata(:,21)/3600)+1.25*(Smoke_N)./(Newdata(:,21)/3600)+2.86*((Smoke_s)*22.4)./(64*Newdata(:,21)/3600)+0.8*((Smoke_W)*22.4)./(18*Newdata(:,21)/3600)+1.429*((Smoke_o2)*22.4)./(32*Newdata(:,21)/3600);
%Ca_pure=0.94;%石灰石纯度
%Cpg=(7.4e-4*(Newdata(:,40)+273)+0.641).*((Smoke_CO2)./(Newdata(:,21)/3600))+(7.4e-4*(Newdata(:,40)+273)+0.641).*((Smoke_s)*22.4./(64*Newdata(:,21)/3600))+(2.5e-4*(Newdata(:,40)+273)+0.845).*((Smoke_o2)*22.4./(32*Newdata(:,21)/3600))+(2.5e-4*(Newdata(:,40)+273)+0.845).*((Smoke_N)./(Newdata(:,21)/3600))+(-5e-4*(Newdata(:,40)+273)+14.6).*((Smoke_W)*22.4./(18*Newdata(:,21)/3600));%烟气比热容kj/kg.K

%烟气比热容
Cpg1=(1.6+1.055e-3*T1-0.696e-6*T1+3.2e-10*T1-0.86e-13*T1).*((Smoke_CO2)./(Newdata(:,21)/3600))...
    +(1.6+1.055e-3*T1-0.696e-6*T1+3.2e-10*T1-0.86e-13*T1).*((Smoke_s)*22.4./(64*Newdata(:,21)/3600))...
    +(1.295-2.25e-5*T1+2.622e-7*T1-1.998e-10*T1+0.646e-13*T1).*((Smoke_o2)*22.4./(32*Newdata(:,21)/3600))...
    +(1.295-2.25e-5*T1+2.622e-7*T1-1.998e-10*T1+0.646e-13*T1).*((Smoke_N)./(Newdata(:,21)/3600))...
    +(1.494+1.038e-4*T1+2.395e-7*T1-1.482e-10*T1+3.833e-14*T1).*((Smoke_W)*22.4./(18*Newdata(:,21)/3600));%烟气比热容 kj/m3.℃ 在炉膛温度及一个大气压下
Cpg2=(1.6+1.055e-3*T2-0.696e-6*T2+3.2e-10*T2-0.86e-13*T2).*((Smoke_CO2)./(Newdata(:,21)/3600))...
    +(1.6+1.055e-3*T2-0.696e-6*T2+3.2e-10*T2-0.86e-13*T2).*((Smoke_s)*22.4./(64*Newdata(:,21)/3600))...
    +(1.295-2.25e-5*T2+2.622e-7*T2-1.998e-10*T2+0.646e-13*T2).*((Smoke_o2)*22.4./(32*Newdata(:,21)/3600))...
    +(1.295-2.25e-5*T2+2.622e-7*T2-1.998e-10*T2+0.646e-13*T2).*((Smoke_N)./(Newdata(:,21)/3600))...
    +(1.494+1.038e-4*T2+2.395e-7*T2-1.482e-10*T2+3.833e-14*T2).*((Smoke_W)*22.4./(18*Newdata(:,21)/3600));%取标况为1am,0℃
% Acon=M_c.*41.79/100;%燃料含灰量kg/s
% Ash_ca=Newdata(:,6).*(1-Ca_pure)*1000/3600;%石灰石杂质 kg/s xx
% CaSO4=(Scon-Smoke_s*32/64)*136/32;%硫酸钙质量kg/s
% Un_CaO=Newdata(:,6).*Ca_pure*0.9*1000/3600-(Scon-Smoke_s*32/64)*56/32;%未反应的CaO kg/s
% Un_ca=(Newdata(:,6).*Ca_pure*0.1*1000)/3600;%未分解的CaCO3 kg/s
% TotalAsh=Acon+Ash_ca+CaSO4+Un_CaO+Un_ca;%总灰量kg/s
% %Unit_ash=TotalAsh./((Newdata(:,5)+Newdata(:,6)+Newdata(:,7))*1000/3600);%单位当量燃料灰分 无量纲
% bottom_ash=TotalAsh-fly_ash;%底灰量kg/s
% Qr=M_c*15880;%单位KJ/s 燃料初始输入热量
% bottom_C=(0.03*Qr./33727-fly_ash*0.0221)./bottom_ash;%底灰含碳量
% Ca_S=mean(((Newdata(577:716,6)*1000/3600).*Ca_pure/100)./(M_c(577:716)*0.372/(32*100)));%钙硫摩尔比
% Q_ca=(0.372*M_c.*(152*desulfu-57.19*2*0.9))/100;%(0.372*M_c.*(501.83*desulfu/32-178.98*Ca_S*0.9/32))/100;石灰石炉内放热KJ/s取自经验公式

P_unit=Newdata(:,67)*1000;%机组功率kW
Kw=mean(P_unit(332:471)./(0.4*(T2(332:471)-Newdata(332:471,45))));%水冷壁换热系数
Q_solid=(Newdata(:,6)*0.59+(Newdata(:,7))*0.92)*1000.*Newdata(:,22)/3600;%输入物料物理显热
%Q_air=Newdata(:,2)*1.29*4.66.*Newdata(:,23)/3600+(1.29*Newdata(:,49).*Newdata(:,47)+1.29*Newdata(:,48).*Newdata(:,46))*4.66/3600;%%一次风和下二次风入炉初始热量
Q_air=(Newdata(:,2)/3600)*371+(Newdata(:,48)/3600)*313.6+(Newdata(:,49)/3600)*362;
Q1=R_c1*15880;%碳燃烧放热
Q_lhf=(Newdata(:,24)*684+Newdata(:,25)*552+Newdata(:,26)*(593)...
+Newdata(:,27)*652+Newdata(:,36)*1135+Newdata(:,37)*1282+Newdata(:,38)*1225+Newdata(:,39)*1192)/3600;
Q_H2O=(Newdata(:,70)*0.312*1000/3600).*(1.77*(273+T1)-4.183*(273+Newdata(:,22)));
%Q_lhf=(Newdata(:,24).*Newdata(:,28)+Newdata(:,25).*Newdata(:,29)+Newdata(:,26).*Newdata(:,30)...
%     +Newdata(:,27).*Newdata(:,31)+Newdata(:,32).*Newdata(:,36)+Newdata(:,33).*Newdata(:,37)...
%     +Newdata(:,34).*Newdata(:,38)+Newdata(:,35).*Newdata(:,39))*1.29*4.66/3600;%流化风带入热量 
Q2=Fu1*0.8.*T1+Fg_m.*Cpg1.*T1;%Fu1和烟气带出热量
%Q_sec=(1.29*Newdata(:,51).*Newdata(:,46)+1.29*Newdata(:,50).*Newdata(:,47))*4.66/3600;%%上二次风带入热量
Q_sec=(Newdata(:,51)*314+Newdata(:,50)*362)/3600;
Q3=R_c2*15880;%稀相区碳燃烧放热
Q4=Fu2.*T2*0.8+Fg.*Cpg2.*T2/3600;%稀相区扬析量和烟气量
fly_ash=Newdata(:,21)*0.17.*Row_gas/(3600*2.21);%飞灰量kg/s
Q5=R_c3*15880;

%外置床
T_hlf1=Newdata(:,32);%回料阀左1温度
T_hlf2=Newdata(:,33);%回料阀左2温度
T_hlf3=Newdata(:,34);%回料阀右1温度
T_hlf4=Newdata(:,35);%回料阀右2温度
T_hlf=(T_hlf1+T_hlf2+T_hlf3+T_hlf4)/4;
h_1=165;%235.98;%kj/kg中过1焓差
h_2=284;%311.68;%中过2焓差
h_3=83;%156.53;%低过焓差
h_4=227.2;%242.26;%高再焓差
T_fl1=Newdata(:,54);%分离器左1温度
T_fl2=Newdata(:,55);%分离器左2温度
T_fl3=Newdata(:,56);%分离器右1温度
T_fl4=Newdata(:,57);%分离器右2温度
T_wch1=Newdata(:,28);%外置床左1出口温度
T_wch2=Newdata(:,30);%外置床左2出口温度
T_wch3=Newdata(:,29);%外置床右1出口温度
T_wch4=Newdata(:,31);%外置床右2出口温度
T_wch=(T_wch1+T_wch2+T_wch3+T_wch4)/4;

Fcold1=((Newdata(:,53)-Newdata(:,59)-Newdata(:,58))*1000*h_1/3600+(Newdata(:,24)/3600)*0.46*1.005.*(T_wch1-Newdata(:,61)))./(T_fl1-T_wch1);%中过1冷灰量
Fcold2=((Newdata(:,53)-Newdata(:,58))*1000*h_2/3600+(Newdata(:,25)/3600)*0.51*1.005.*(T_wch3-Newdata(:,63)))./(T_fl3-T_wch3);%中过2冷灰量
Fcold3=((Newdata(:,53)-Newdata(:,59)-Newdata(:,58)-Newdata(:,60))*1000*h_3/3600+(Newdata(:,26)/3600)*0.49*1.005.*(T_wch2-Newdata(:,62)))./(T_fl2-T_wch2);%外置床B冷灰
Fcold4=((Newdata(:,53))*1000*h_4/3600+(Newdata(:,27)/3600)*0.46*1.005.*(T_wch4-Newdata(:,64)))./(T_fl4-T_wch4);%外置床D冷灰量kg
F_cold=Fcold1+Fcold2+Fcold3+Fcold4;
% %Fcold1=(889502*h_1/3600+(Newdata(:,24)/3600)*1.29*4.66.*(T_wch1-Newdata(:,61)))./(0.8*(T_fl1-T_wch1));%%中过1冷灰量
% fcold1=(889502*h_1/3600+(Newdata(:,24)/3600)*673.5)./(0.8*(T_fl1-T_wch1));
% %Fcold2=(929465*h_2/3600+(Newdata(:,25)/3600)*1.29*4.66.*(T_wch3-Newdata(:,63)))./(0.8*(T_fl3-T_wch3));%%中过2冷灰量
% fcold2=(929465*h_2/3600+(Newdata(:,25)/3600)*539)./(0.8*(T_fl3-T_wch3));
% %Fcold3=(887080*h_3/3600+(Newdata(:,26)/3600)*1.29*4.66.*(T_wch2-Newdata(:,62)))./(0.8*(T_fl2-T_wch2));%%%外置床B冷灰
% fcold3=(887080*h_3/3600+(Newdata(:,26)/3600)*583)./(0.8*(T_fl2-T_wch2));
% %Fcold4=(807054*h_4/3600+(Newdata(:,27)/3600)*1.29*4.66.*(T_wch4-Newdata(:,64)))./(0.8*(T_fl4-T_wch4));%%外置床D冷灰量kg
% fcold4=(807054*h_4/3600+(Newdata(:,27)/3600)*640)./(0.8*(T_fl4-T_wch4));
% F_cold=fcold1+fcold2+fcold3+fcold4;%Fcold1+Fcold2+Fcold3+Fcold4;%冷灰量
Q_cold=0.8*(Fcold1.*T_wch1+Fcold2.*T_wch3+Fcold3.*T_wch2+Fcold4.*T_wch4);%冷灰带入热量
Q_hot=0.8*(Fcc-F_cold).*(T_hlf);%热灰带入热量
Q_cc=Q_cold+Q_hot;%循环物料带入热量

%数据验证
f1=Q_solid(332:471)+Q_air(332:471)+Q1(332:471)+Q_lhf(332:471)+Q_cc(332:471)+Q_hff(332:471)-Q_H2O(332:471);
%f2=Fu1(2000:3000)*0.8;
f21=Q2(332:471);%((273+T1(2000:3000))/273).*Fg_m(2000:3000).*Cpg1(2000:3000);
%f22=Q_cold(2300:3000);
%f23=T_wch(2300:3000);
f3=Q_sec(332:471)+Q3(332:471);
%f4=Fu2(2000:3000)*0.8;
f41=Q4(332:471);%((273+T2(2000:3000))/273).*Newdata(2000:3000,21).*Cpg2(2000:3000)/3600;
f5=Kw;
f6=Newdata(332:471,45);
f7=M1(332:471)*0.8;
f8=M2(332:471)*0.8;
f9=Newdata(332:471,21).*Cpg2(332:471)/3600;
f91=Q5(332:471);
f92=Fcc(332:471)*0.8;
f93=fly_ash(332:471)*0.8;
f11=Fu2(332:471)*0.8;
out=zeros(140,3);
out(1,:)=[904,854,906];
for k=1:1:139
  out(k+1,1)=(f1(k)-f21(k))/f7(k)+out(k,1);
  out(k+1,2)=(f21(k)+f3(k)-f41(k)-f5*(out(k,2)-f6(k)))/f8(k)+out(k,2);
  out(k+1,3)=(f41(k)+f91(k)-f92(k)*out(k,3)-f93(k)*out(k,3)-f9(k)*out(k,3))/f11(k)+out(k,3);
end
out
x1=[332:1:471];
L=out(:,1);
% q=R_c1(2300:3000);
T=Newdata(332:471,9);
figure(2)
plot(x1,L)
hold on
plot(x1,T,'k')
% figure(3)
% plot(x1,T2(2300:3000),'r')
% figure(3)
% plot(x1,q,'b')
% figure(4)
% plot(x1,T3(2300:3000))




