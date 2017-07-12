Newdata=xlsread('1','sheet1');   
Sb=117.8;%密相区通风截面积
Sv=229.15;% m2
h1=4.78;%密相区压力测点之间距离 m
h2=6.8;%稀相区测点至炉膛高度距离 m
Fg_m=(Newdata(:,2)+Newdata(:,8)+Newdata(:,3))/3600;%密相区风量Nm3/s  
Fg_x=(Newdata(:,1))/3600;%稀相区上二次风 Nm3/s
Fg=Newdata(:,11)/3600;%总风量 Nm3/s
Tg1=Newdata(:,23);%空预器出口一次风温度
T1=Newdata(:,9);%密相区床温
T2=Newdata(:,41);%稀相区床温
V1=Fg_m.*(273+T1)/(273*Sb);%密相区实际流化速度m/s
V2=Fg.*(273+T2)/(273*Sv);%稀相区实际流化速度m/s
oxy_air=Newdata(:,42)/100;%空预器出口含氧量“%”
dc1=0.002;row_c1=2000;%密相区颗粒直径以及颗粒密度
dc2=0.0002;row_c2=1800;
Vb=1796.1;
Vs=4033.1;%各相区体积m3
delta_p1=(Newdata(:,43)-Newdata(:,44))*1000;%密相区床压差，根据测点图确定
delta_p2=Newdata(:,15);%炉膛上部差压
M1=delta_p1*Vb/(9.8*h1);
M2=delta_p2*Vs/(9.8*h2);%各相区床料质量

row_g=1.293;%空气密度在0℃，标准大气压下
row_gm1=row_g*(273./(273+T1)).*((Newdata(:,43)*1000+101325)/101325);%在密相区的温度压力下，空气密度
row_gm2=row_g*(273./(273+T2)).*((150+101325)/101325);

Ks1=0.513*(T1+273).*exp(-9160./(273+T1));%燃烧速率Rc1常数
Ks2=0.513*(T2+273).*exp(-9160./(T2+273));%燃烧速率Rc2常数
p=2500*exp(-5.19e-4./(8.319*(273+120+T2)));
phai=(2*p-(p/0.095)*(100*dc2-0.005))./(p+2);
Dg1=(3.13e-4*((273+T1)/1500).^1.75).*(101./(Newdata(:,43)+101));%氧气扩散系数 m2/s 考虑为炉膛下部压力
Dg2=(3.13e-4*((273+T2)/1500).^1.75).*(101325/(150+101325));%炉膛上部压力测点待定，影响计算
U_d1=17.9e-6*(((273+T1)/298).^1.5).*(408.4./(383.4+T1));%1.544e-4*row_gm1;%动力粘度
U_d2=17.9e-6*(((273+T2)/298).^1.5).*(408.4./(383.4+T2));
v1=U_d1.*row_gm1;%运动粘度
v2=U_d2.*row_gm2;
Sc1=U_d1./(row_gm1.*Dg1);%Schmidt数
Sc2=U_d2./(row_gm2.*Dg2);
Ar1=(row_gm1*(dc1)^3.*(row_c1-row_gm1)*9.8)./(row_gm1*1.544e-4).^2;%阿基米德数
Ar2=(row_gm2*(dc2)^3.*(row_c2-row_gm2)*9.8)./(row_gm2*1.421e-4).^2;
e1=1-M1./(Vb*(row_c1-row_gm1));%密相区空隙率
e2=1-M2./(Vs*(row_c2-row_gm2));
Re1=(Ar1.*e1.^4.75)./(18+0.61*(Ar1.*e1.^4.74).^0.5);%雷诺数
Re2=(Ar2.*e2.^4.75)./(18+0.61*(Ar2.*e2.^4.74).^0.5);
Sh1=2+0.6*(Re1.^0.5).*(Sc1.^0.33);%舍伍德数
Sh2=2+0.6*(Re2.^0.5).*(Sc2.^0.33);
Kd1=Sh1.*Dg1/dc1;%扩散反应速率
Kd2=Sh2.*Dg2/dc2;
kc1=1./(1./Ks1+1./Kd1);%燃烧速率Rc1常数 m/s
kc2=1./(1./Ks2+1./(phai.*Kd2));

C_o2=mean(oxy_air/22.4);%稀相区氧气摩尔浓度kmol/m3
C_o1=0.09;%(Newdata(:,1)+Newdata(:,2))/3600-R_c1*22.4/12;%密相区氧气浓度kg/m3


% R_c1=72*Kc1*C_o1*x/(0.002*1800);
% R_c2=72*Kc2*C_o2*x2/(0.0002*1800);
% Fu1=y*M1.*(V1-4.5);
% Fu2=y2*M2.*(V2-4.5);
M=(Newdata(:,5)+Newdata(:,6)+Newdata(:,7))*1000/3600;%单位时间床料总量kg/s
M_c=(Newdata(:,5)+Newdata(:,7))*1000/3600;%单位时间内给煤量 kg/s
F_hff=(Newdata(:,5)+Newdata(:,7))*173.25/3600;%取自空干基混煤挥发分平均值 kg/s
N_fl=Newdata(:,17);%分离效率
a=30*M_c;
Vmf1=0.294*(0.002^0.584./(v1).^0.056).*((row_c1-row_gm1)./row_gm1).^0.528;%临界流化速度 row_gm1
Vmf2=0.294*(0.0002^0.584./(v2).^0.056).*((row_c2-row_gm2)./row_gm2).^0.528;
C_bed=(Newdata(:,5)+Newdata(:,7))*1000*0.4103/3600;%燃料量中碳含量 kg/s
Fg_gas=Newdata(:,21).*Newdata(:,20)/(100*3600);%烟气含氧量 Nm3/s
Y=zeros(481,6);

%while C_o11<0 %|| abs(C_o11-C_o1)>0.01
    %C_o1=(C_o1+C_o11)/2
%     if C_o11>C_o1
%     C_o1=C_o1-(C_o11-C_o1)/2;
%     else
%       C_o1=C_o1+(C_o1-C_o11)/2;  
%     end
for i=332:1:812
   x0=ones(1,6);
   a1=M(i)-F_hff(i);
   a2=N_fl(i)*M2(i)*(V2(i)-Vmf2(i));
   a3=72*kc1(i)/(dc1*row_c1);
   a4=M1(i)*(V1(i)-Vmf1(i));
   b1=72*kc2(i)*C_o2/(dc2*row_c2);
   b2=M2(i)*(V2(i)-Vmf2(i));
   c1=C_bed(i);
   c2=N_fl(i)*(V2(i)-Vmf2(i));
   c3=V1(i)-Vmf1(i);
   d1=V2(i)-Vmf2(i);
   d2=0.21*Fg_m(i)/22.4;
   d3=0.21*Fg_x(i)/22.4;
   d4=Fg_gas(i)/22.4;
   x=fsolve(@(x) myfun(x,a1,a2,a3,a4,b1,b2,c1,c2,c3,d1,d2,d3,d4),x0);
   Y(i-331,:)=x;
end
Y
y=mean(Y);
m_c1=y(1);
m_c2=y(2);
k1=y(3);
k2=y(4);
C_o1=y(5);
Fg_out1=y(6);
R_c1=72*C_o1*kc1*m_c1/(dc1*row_c1);%燃烧速率kg/s
R_c2=72*C_o2*kc2.*m_c2/(dc2*row_c2);
Fu1=k1*M1.*(V1-Vmf1);%密相区进入稀相区的物料量kg/s
Fu2=k2*M2.*(V2-Vmf2);
Fcc=N_fl.*Fu2;%循环物料量
%Uv_o2=0.21*Newdata(:,11)/3600-((R_c1+R_c2)*22.4/12)-(Newdata(:,21).*Newdata(:,20)/3600);%挥发份消耗氧量m3/s
% sum(Newdata(332:812,3)*0.21/3600)
% %sum((Newdata(332:812,1)+Newdata(332:812,2))*0.21/3600)
% sum(R_c1*22.4/12)
% F_out1=C_o2*Vs+1.429*(sum((Newdata(332:812,20).*Newdata(332:812,21))/(100*60))-sum(Newdata(332:812,3)*0.21/60))+sum(R_c2*60*32/12);%由密相区进入稀相区的氧气含量
% C_o11=(1.249*sum((Newdata(332:812,1)+Newdata(332:812,2))*0.21/60)-F_out1-sum(R_c1*60*32/12))/Vb;

Scon=M_c*0.3722/100;%燃料含S量kg/s
Smoke_s=Newdata(:,21).*Newdata(:,18)/(1000000*3600);%烟气含SO2量 kg/s
desulfu=(Scon-Smoke_s*32/64)./Scon;%脱硫效率 矩阵
Smoke_o2=Fg_gas*32/22.4;%氧气单位kg/s
%Smoke_NOx=Newdata(:,21).*Newdata(:,19)/(1000000*3600);%烟气含NOx量kg/s
Smoke_N=M_c*0.7171*22.4/(100*28)+Newdata(:,11)*0.78/3600;%烟气含N2量m3/s
Smoke_W=Newdata(:,11)*7.23/(1000*3600)+M_c*1.38/100+M_c*7.51/100;%烟气含水量kg/s
Smoke_CO2=Newdata(:,21)/3600-Smoke_s*22.4/64-Smoke_o2*22.4/32-Smoke_N-Smoke_W*22.4/18;%烟气CO2含量m3/s
row_gas=1.96*(Smoke_CO2)./(Newdata(:,21)/3600)+1.25*(Smoke_N)./(Newdata(:,21)/3600)+2.86*((Smoke_s)*22.4)./(64*Newdata(:,21)/3600)+0.8*((Smoke_W)*22.4)./(18*Newdata(:,21)/3600)+1.429*((Smoke_o2)*22.4)./(32*Newdata(:,21)/3600);%烟气密度 kg/m3
Ca_pure=0.95;%(Scon-Smoke_s*32/64)*100)./(Newdata(:,6).*desulfu*0.9*32*1000/3600);%石灰石纯度
Cpg=(0.00032*Newdata(:,40)+0.912).*((Smoke_CO2)./(Newdata(:,21)/3600))+(0.00032*Newdata(:,40)+0.912).*((Smoke_s)*22.4./(64*Newdata(:,21)/3600))+(0.00175*Newdata(:,40)+0.91).*((Smoke_o2)*22.4./(32*Newdata(:,21)/3600))+(0.00175*Newdata(:,40)+0.91).*((Smoke_N)./(Newdata(:,21)/3600))+(0.00115*Newdata(:,40)+13.61).*((Smoke_W)*22.4./(18*Newdata(:,21)/3600));%烟气比热容

Acon=M_c.*41.79/100;%燃料含灰量kg/s
Ash_ca=Newdata(:,6).*(1-Ca_pure)*1000/3600;%石灰石杂质 kg/s xx
CaSO4=(Scon-Smoke_s*32/64)*136/32;%硫酸钙质量kg/s
Un_CaO=Newdata(:,6).*Ca_pure*0.9*1000/3600-(Scon-Smoke_s*32/64)*56/32;%未反应的CaO kg/s
Un_ca=(Newdata(:,6).*Ca_pure*0.1*1000)/3600;%未分解的CaCO3 kg/s
TotalAsh=Acon+Ash_ca+CaSO4+Un_CaO+Un_ca;%总灰量kg/s
%Unit_ash=TotalAsh./((Newdata(:,5)+Newdata(:,6)+Newdata(:,7))*1000/3600);%单位当量燃料灰分 无量纲
fly_ash=Newdata(:,21)*0.17.*row_gas/3600*2.21;%飞灰量kg/s
bottom_ash=TotalAsh-fly_ash;%底灰量kg/s
Qr=M_c*15880;%单位KJ/s 燃料初始输入热量
bottom_C=(0.03*Qr./33727-fly_ash*0.0221)./bottom_ash;%底灰含碳量

T_hlf1=Newdata(:,32);%回料阀左1温度
T_hlf2=Newdata(:,33);%回料阀左2温度
T_hlf3=Newdata(:,34);%回料阀右1温度
T_hlf4=Newdata(:,35);%回料阀右2温度
Ca_S=((Newdata(:,6)*1000/3600)*Ca_pure/100)./(M_c*0.372/(32*100));%钙硫摩尔比
Q_ca=(0.372*M_c.*(152*desulfu-57.19*Ca_S*0.9))/100;%石灰石炉内放热KJ/s 取自经验公式
Q_solid=(Newdata(:,6)*0.59+(Newdata(:,5)+Newdata(:,7))*0.92)*1000.*Newdata(:,22)/3600;%输入物料物理显热
Q_air=Newdata(:,2)*0.64*1.005.*Newdata(:,23)/3600+(0.65*Newdata(:,49).*Newdata(:,47)+0.69*Newdata(:,48).*Newdata(:,46))*1.005/3600;%一次风和下二次风入炉初始热量
Q1=R_c1*15880;%碳燃烧放热
Q_lhf=(0.46*Newdata(:,24).*Newdata(:,28)+0.51*Newdata(:,25).*Newdata(:,29)+0.49*Newdata(:,26).*Newdata(:,30)...
    +0.46*Newdata(:,27).*Newdata(:,31)+0.33*Newdata(:,32).*Newdata(:,36)+0.3*Newdata(:,33).*Newdata(:,37)...
    +0.3*Newdata(:,34).*Newdata(:,38)+0.31*Newdata(:,35).*Newdata(:,39))*1.005/3600;%流化风带入热量 气体密度按平均温度计算
Q2=Fu1*0.8.*T1+Fg_m.*Cpg.*T1;%Fu1和烟气带出热量
Q_sec=(0.69*Newdata(:,51).*Newdata(:,46)+0.65*Newdata(:,50).*Newdata(:,47))*1.005/3600;%上二次风带入热量
Q3=R_c2*15880;%稀相区碳燃烧放热
Q4=Fu2.*T2*0.8+row_gas.*Newdata(:,21).*Cpg.*T2/3600;%稀相区扬析量和烟气量

Kw=mean((Q2(332:812)+Q_sec(332:812)+Q3(332:812)-Q4(332:812))./(1056.6*(T2(332:812)-Newdata(332:812,45))));%水冷壁换热系数
h_1=165;%kj/kg中过1焓差
h_2=284;%中过2焓差
h_3=83;%低过焓差
h_4=227.8;%高再焓差
T_fl1=Newdata(:,54);%分离器左1温度
T_fl2=Newdata(:,55);%分离器左2温度
T_fl3=Newdata(:,56);%分离器右1温度
T_fl4=Newdata(:,57);%分离器右2温度
T_wch1=Newdata(:,28);%外置床左1出口温度
T_wch2=Newdata(:,30);%外置床左2出口温度
T_wch3=Newdata(:,29);%外置床右1出口温度
T_wch4=Newdata(:,31);%外置床右2出口温度
Fcold1=((Newdata(:,53)-Newdata(:,59)-Newdata(:,58))*1000*h_1/3600+(Newdata(:,24)/3600)*0.46*1.005.*(T_wch1-Newdata(:,61)))./(T_fl1-T_wch1);%中过1冷灰量
Fcold2=((Newdata(:,53)-Newdata(:,58))*1000*h_2/3600+(Newdata(:,25)/3600)*0.51*1.005.*(T_wch3-Newdata(:,63)))./(T_fl3-T_wch3);%中过2冷灰量
Fcold3=((Newdata(:,53)-Newdata(:,59)-Newdata(:,58)-Newdata(:,60))*1000*h_3/3600+(Newdata(:,26)/3600)*0.49*1.005.*(T_wch2-Newdata(:,62)))./(T_fl2-T_wch2);%外置床B冷灰
Fcold4=((Newdata(:,53))*1000*h_4/3600+(Newdata(:,27)/3600)*0.46*1.005.*(T_wch4-Newdata(:,64)))./(T_fl4-T_wch4);%外置床D冷灰量kg
F_cold=Fcold1+Fcold2+Fcold3+Fcold4;%冷灰量
Q_cold=0.8*(Fcold1.*T_wch1+Fcold2.*T_wch3+Fcold3.*T_wch2+Fcold4.*T_wch4);%冷灰带入热量
Q_hot=0.8*(Fcc-F_cold).*((T_hlf1+T_hlf2+T_hlf3+T_hlf4)/4);%热灰带入热量
Q_cc=Q_cold+Q_hot;%循环物料带入热量
%数据验证

f1=Q_solid(2000:3000)+Q_air(2000:3000)+Q1(2000:3000)+Q_lhf(2000:3000)+Q_ca(2000:3000)+Q_cc(2000:3000);
f2=Fu1(2000:3000)*0.8+Fg_m(2000:3000).*Cpg(2000:3000);
f3=Q_sec(2000:3000)+Q3(2000:3000);
f4=Fu2(2000:3000)*0.8+row_gas(2000:3000).*Newdata(2000:3000,21).*Cpg(2000:3000)/3600;
f5=Kw*1056.6;
f6=Newdata(2000:3000,45);
f7=M1(2000:3000)*0.8;
f8=M2(2000:3000)*0.8;
out=zeros(1000,2);
out(1,:)=[900,823];
for k=1:1:1000
  out(k+1,1)=(f1(k)+f7(k)*out(k,1))/(f7(k)+f2(k));
  out(k+1,2)=(f2(k)*out(k+1,1)+f3(k)+f5*f6(k)+f8(k)*out(k,2))/(f4(k)+f5+f8(k));
end
out





