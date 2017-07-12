Newdata=xlsread('1','sheet1');   
Sb=117.8;%������ͨ������
Sv=229.15;% m2
h1=4.78;%������ѹ�����֮����� m
h2=6.8;%ϡ���������¯�Ÿ߶Ⱦ��� m
Fg_m=(Newdata(:,2)+Newdata(:,8)+Newdata(:,3))/3600;%����������Nm3/s  
Fg_x=(Newdata(:,1))/3600;%ϡ�����϶��η� Nm3/s
Fg=Newdata(:,11)/3600;%�ܷ��� Nm3/s
Tg1=Newdata(:,23);%��Ԥ������һ�η��¶�
T1=Newdata(:,9);%����������
T2=Newdata(:,41);%ϡ��������
V1=Fg_m.*(273+T1)/(273*Sb);%������ʵ�������ٶ�m/s
V2=Fg.*(273+T2)/(273*Sv);%ϡ����ʵ�������ٶ�m/s
oxy_air=Newdata(:,42)/100;%��Ԥ�����ں�������%��
dc1=0.002;row_c1=2000;%����������ֱ���Լ������ܶ�
dc2=0.0002;row_c2=1800;
Vb=1796.1;
Vs=4033.1;%���������m3
delta_p1=(Newdata(:,43)-Newdata(:,44))*1000;%��������ѹ����ݲ��ͼȷ��
delta_p2=Newdata(:,15);%¯���ϲ���ѹ
M1=delta_p1*Vb/(9.8*h1);
M2=delta_p2*Vs/(9.8*h2);%��������������

row_g=1.293;%�����ܶ���0�棬��׼����ѹ��
row_gm1=row_g*(273./(273+T1)).*((Newdata(:,43)*1000+101325)/101325);%�����������¶�ѹ���£������ܶ�
row_gm2=row_g*(273./(273+T2)).*((150+101325)/101325);

Ks1=0.513*(T1+273).*exp(-9160./(273+T1));%ȼ������Rc1����
Ks2=0.513*(T2+273).*exp(-9160./(T2+273));%ȼ������Rc2����
p=2500*exp(-5.19e-4./(8.319*(273+120+T2)));
phai=(2*p-(p/0.095)*(100*dc2-0.005))./(p+2);
Dg1=(3.13e-4*((273+T1)/1500).^1.75).*(101./(Newdata(:,43)+101));%������ɢϵ�� m2/s ����Ϊ¯���²�ѹ��
Dg2=(3.13e-4*((273+T2)/1500).^1.75).*(101325/(150+101325));%¯���ϲ�ѹ����������Ӱ�����
U_d1=17.9e-6*(((273+T1)/298).^1.5).*(408.4./(383.4+T1));%1.544e-4*row_gm1;%����ճ��
U_d2=17.9e-6*(((273+T2)/298).^1.5).*(408.4./(383.4+T2));
v1=U_d1.*row_gm1;%�˶�ճ��
v2=U_d2.*row_gm2;
Sc1=U_d1./(row_gm1.*Dg1);%Schmidt��
Sc2=U_d2./(row_gm2.*Dg2);
Ar1=(row_gm1*(dc1)^3.*(row_c1-row_gm1)*9.8)./(row_gm1*1.544e-4).^2;%�����׵���
Ar2=(row_gm2*(dc2)^3.*(row_c2-row_gm2)*9.8)./(row_gm2*1.421e-4).^2;
e1=1-M1./(Vb*(row_c1-row_gm1));%��������϶��
e2=1-M2./(Vs*(row_c2-row_gm2));
Re1=(Ar1.*e1.^4.75)./(18+0.61*(Ar1.*e1.^4.74).^0.5);%��ŵ��
Re2=(Ar2.*e2.^4.75)./(18+0.61*(Ar2.*e2.^4.74).^0.5);
Sh1=2+0.6*(Re1.^0.5).*(Sc1.^0.33);%�������
Sh2=2+0.6*(Re2.^0.5).*(Sc2.^0.33);
Kd1=Sh1.*Dg1/dc1;%��ɢ��Ӧ����
Kd2=Sh2.*Dg2/dc2;
kc1=1./(1./Ks1+1./Kd1);%ȼ������Rc1���� m/s
kc2=1./(1./Ks2+1./(phai.*Kd2));

C_o2=mean(oxy_air/22.4);%ϡ��������Ħ��Ũ��kmol/m3
C_o1=0.09;%(Newdata(:,1)+Newdata(:,2))/3600-R_c1*22.4/12;%����������Ũ��kg/m3


% R_c1=72*Kc1*C_o1*x/(0.002*1800);
% R_c2=72*Kc2*C_o2*x2/(0.0002*1800);
% Fu1=y*M1.*(V1-4.5);
% Fu2=y2*M2.*(V2-4.5);
M=(Newdata(:,5)+Newdata(:,6)+Newdata(:,7))*1000/3600;%��λʱ�䴲������kg/s
M_c=(Newdata(:,5)+Newdata(:,7))*1000/3600;%��λʱ���ڸ�ú�� kg/s
F_hff=(Newdata(:,5)+Newdata(:,7))*173.25/3600;%ȡ�Կոɻ���ú�ӷ���ƽ��ֵ kg/s
N_fl=Newdata(:,17);%����Ч��
a=30*M_c;
Vmf1=0.294*(0.002^0.584./(v1).^0.056).*((row_c1-row_gm1)./row_gm1).^0.528;%�ٽ������ٶ� row_gm1
Vmf2=0.294*(0.0002^0.584./(v2).^0.056).*((row_c2-row_gm2)./row_gm2).^0.528;
C_bed=(Newdata(:,5)+Newdata(:,7))*1000*0.4103/3600;%ȼ������̼���� kg/s
Fg_gas=Newdata(:,21).*Newdata(:,20)/(100*3600);%���������� Nm3/s
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
R_c1=72*C_o1*kc1*m_c1/(dc1*row_c1);%ȼ������kg/s
R_c2=72*C_o2*kc2.*m_c2/(dc2*row_c2);
Fu1=k1*M1.*(V1-Vmf1);%����������ϡ������������kg/s
Fu2=k2*M2.*(V2-Vmf2);
Fcc=N_fl.*Fu2;%ѭ��������
%Uv_o2=0.21*Newdata(:,11)/3600-((R_c1+R_c2)*22.4/12)-(Newdata(:,21).*Newdata(:,20)/3600);%�ӷ�����������m3/s
% sum(Newdata(332:812,3)*0.21/3600)
% %sum((Newdata(332:812,1)+Newdata(332:812,2))*0.21/3600)
% sum(R_c1*22.4/12)
% F_out1=C_o2*Vs+1.429*(sum((Newdata(332:812,20).*Newdata(332:812,21))/(100*60))-sum(Newdata(332:812,3)*0.21/60))+sum(R_c2*60*32/12);%������������ϡ��������������
% C_o11=(1.249*sum((Newdata(332:812,1)+Newdata(332:812,2))*0.21/60)-F_out1-sum(R_c1*60*32/12))/Vb;

Scon=M_c*0.3722/100;%ȼ�Ϻ�S��kg/s
Smoke_s=Newdata(:,21).*Newdata(:,18)/(1000000*3600);%������SO2�� kg/s
desulfu=(Scon-Smoke_s*32/64)./Scon;%����Ч�� ����
Smoke_o2=Fg_gas*32/22.4;%������λkg/s
%Smoke_NOx=Newdata(:,21).*Newdata(:,19)/(1000000*3600);%������NOx��kg/s
Smoke_N=M_c*0.7171*22.4/(100*28)+Newdata(:,11)*0.78/3600;%������N2��m3/s
Smoke_W=Newdata(:,11)*7.23/(1000*3600)+M_c*1.38/100+M_c*7.51/100;%������ˮ��kg/s
Smoke_CO2=Newdata(:,21)/3600-Smoke_s*22.4/64-Smoke_o2*22.4/32-Smoke_N-Smoke_W*22.4/18;%����CO2����m3/s
row_gas=1.96*(Smoke_CO2)./(Newdata(:,21)/3600)+1.25*(Smoke_N)./(Newdata(:,21)/3600)+2.86*((Smoke_s)*22.4)./(64*Newdata(:,21)/3600)+0.8*((Smoke_W)*22.4)./(18*Newdata(:,21)/3600)+1.429*((Smoke_o2)*22.4)./(32*Newdata(:,21)/3600);%�����ܶ� kg/m3
Ca_pure=0.95;%(Scon-Smoke_s*32/64)*100)./(Newdata(:,6).*desulfu*0.9*32*1000/3600);%ʯ��ʯ����
Cpg=(0.00032*Newdata(:,40)+0.912).*((Smoke_CO2)./(Newdata(:,21)/3600))+(0.00032*Newdata(:,40)+0.912).*((Smoke_s)*22.4./(64*Newdata(:,21)/3600))+(0.00175*Newdata(:,40)+0.91).*((Smoke_o2)*22.4./(32*Newdata(:,21)/3600))+(0.00175*Newdata(:,40)+0.91).*((Smoke_N)./(Newdata(:,21)/3600))+(0.00115*Newdata(:,40)+13.61).*((Smoke_W)*22.4./(18*Newdata(:,21)/3600));%����������

Acon=M_c.*41.79/100;%ȼ�Ϻ�����kg/s
Ash_ca=Newdata(:,6).*(1-Ca_pure)*1000/3600;%ʯ��ʯ���� kg/s xx
CaSO4=(Scon-Smoke_s*32/64)*136/32;%���������kg/s
Un_CaO=Newdata(:,6).*Ca_pure*0.9*1000/3600-(Scon-Smoke_s*32/64)*56/32;%δ��Ӧ��CaO kg/s
Un_ca=(Newdata(:,6).*Ca_pure*0.1*1000)/3600;%δ�ֽ��CaCO3 kg/s
TotalAsh=Acon+Ash_ca+CaSO4+Un_CaO+Un_ca;%�ܻ���kg/s
%Unit_ash=TotalAsh./((Newdata(:,5)+Newdata(:,6)+Newdata(:,7))*1000/3600);%��λ����ȼ�ϻҷ� ������
fly_ash=Newdata(:,21)*0.17.*row_gas/3600*2.21;%�ɻ���kg/s
bottom_ash=TotalAsh-fly_ash;%�׻���kg/s
Qr=M_c*15880;%��λKJ/s ȼ�ϳ�ʼ��������
bottom_C=(0.03*Qr./33727-fly_ash*0.0221)./bottom_ash;%�׻Һ�̼��

T_hlf1=Newdata(:,32);%���Ϸ���1�¶�
T_hlf2=Newdata(:,33);%���Ϸ���2�¶�
T_hlf3=Newdata(:,34);%���Ϸ���1�¶�
T_hlf4=Newdata(:,35);%���Ϸ���2�¶�
Ca_S=((Newdata(:,6)*1000/3600)*Ca_pure/100)./(M_c*0.372/(32*100));%����Ħ����
Q_ca=(0.372*M_c.*(152*desulfu-57.19*Ca_S*0.9))/100;%ʯ��ʯ¯�ڷ���KJ/s ȡ�Ծ��鹫ʽ
Q_solid=(Newdata(:,6)*0.59+(Newdata(:,5)+Newdata(:,7))*0.92)*1000.*Newdata(:,22)/3600;%����������������
Q_air=Newdata(:,2)*0.64*1.005.*Newdata(:,23)/3600+(0.65*Newdata(:,49).*Newdata(:,47)+0.69*Newdata(:,48).*Newdata(:,46))*1.005/3600;%һ�η���¶��η���¯��ʼ����
Q1=R_c1*15880;%̼ȼ�շ���
Q_lhf=(0.46*Newdata(:,24).*Newdata(:,28)+0.51*Newdata(:,25).*Newdata(:,29)+0.49*Newdata(:,26).*Newdata(:,30)...
    +0.46*Newdata(:,27).*Newdata(:,31)+0.33*Newdata(:,32).*Newdata(:,36)+0.3*Newdata(:,33).*Newdata(:,37)...
    +0.3*Newdata(:,34).*Newdata(:,38)+0.31*Newdata(:,35).*Newdata(:,39))*1.005/3600;%������������� �����ܶȰ�ƽ���¶ȼ���
Q2=Fu1*0.8.*T1+Fg_m.*Cpg.*T1;%Fu1��������������
Q_sec=(0.69*Newdata(:,51).*Newdata(:,46)+0.65*Newdata(:,50).*Newdata(:,47))*1.005/3600;%�϶��η��������
Q3=R_c2*15880;%ϡ����̼ȼ�շ���
Q4=Fu2.*T2*0.8+row_gas.*Newdata(:,21).*Cpg.*T2/3600;%ϡ������������������

Kw=mean((Q2(332:812)+Q_sec(332:812)+Q3(332:812)-Q4(332:812))./(1056.6*(T2(332:812)-Newdata(332:812,45))));%ˮ��ڻ���ϵ��
h_1=165;%kj/kg�й�1�ʲ�
h_2=284;%�й�2�ʲ�
h_3=83;%�͹��ʲ�
h_4=227.8;%�����ʲ�
T_fl1=Newdata(:,54);%��������1�¶�
T_fl2=Newdata(:,55);%��������2�¶�
T_fl3=Newdata(:,56);%��������1�¶�
T_fl4=Newdata(:,57);%��������2�¶�
T_wch1=Newdata(:,28);%���ô���1�����¶�
T_wch2=Newdata(:,30);%���ô���2�����¶�
T_wch3=Newdata(:,29);%���ô���1�����¶�
T_wch4=Newdata(:,31);%���ô���2�����¶�
Fcold1=((Newdata(:,53)-Newdata(:,59)-Newdata(:,58))*1000*h_1/3600+(Newdata(:,24)/3600)*0.46*1.005.*(T_wch1-Newdata(:,61)))./(T_fl1-T_wch1);%�й�1�����
Fcold2=((Newdata(:,53)-Newdata(:,58))*1000*h_2/3600+(Newdata(:,25)/3600)*0.51*1.005.*(T_wch3-Newdata(:,63)))./(T_fl3-T_wch3);%�й�2�����
Fcold3=((Newdata(:,53)-Newdata(:,59)-Newdata(:,58)-Newdata(:,60))*1000*h_3/3600+(Newdata(:,26)/3600)*0.49*1.005.*(T_wch2-Newdata(:,62)))./(T_fl2-T_wch2);%���ô�B���
Fcold4=((Newdata(:,53))*1000*h_4/3600+(Newdata(:,27)/3600)*0.46*1.005.*(T_wch4-Newdata(:,64)))./(T_fl4-T_wch4);%���ô�D�����kg
F_cold=Fcold1+Fcold2+Fcold3+Fcold4;%�����
Q_cold=0.8*(Fcold1.*T_wch1+Fcold2.*T_wch3+Fcold3.*T_wch2+Fcold4.*T_wch4);%��Ҵ�������
Q_hot=0.8*(Fcc-F_cold).*((T_hlf1+T_hlf2+T_hlf3+T_hlf4)/4);%�ȻҴ�������
Q_cc=Q_cold+Q_hot;%ѭ�����ϴ�������
%������֤

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





