clc;
%% Antoine Parameters
A = 6.93666;
B = 1460.793;
C = 207.777;
%% Critical values for cumene
Tc = 631; %% Kelvin
Pc = 32.09; %% in bar
R = 0.08314;
%% Parameters for Peng Robinson equation of state
a = 0.45724*(R^2)*(Tc^2)/Pc;
b = 0.0778*R*Tc/Pc;
w = 0.3274;
sigma = 2.4142;
eps = -0.4142;

i = 1;
P_array = zeros(20,1);   %% store pressure values
vg = zeros(20,1);        %% store maximum volume
vl = zeros(20,1);        %% store minimum volume
for T = 400:10:600
    
    Psat = 10^(A-B/(T-273.15+C)); %% in mmHg
    Psat = Psat/760; %% in bar
    alpha = (1+(0.37464+1.54226*w-0.26992*w^2)*(1-sqrt(T/Tc)))^2;
    v_array = linspace(0,10,1000);
    pressure = -(R.*T./(v_array-b)-a.*alpha./(v_array.^2+2.*b.*v_array-b^2));
    plot(v_array,pressure);
    xlabel("Volume (in L/mol)");
    ylabel("Pressure (in bar)");
    grid on;
    hold on;
    
    a11 = R*T*(b^2)+Psat*(b^3)-a*alpha*b;
    a21 = a*alpha-2*R*T*b-3*Psat*(b^2);
    a31 = Psat*b-R*T;
    a41 = Psat;
    coeff_mat = [a41 a31 a21 a11];     %% Equation is converted to Pv^3+(Pb-RT)v^2 +v(a*aplha-2bRT-3Pb^2)+(Pb^3-a*aplha*b+RTb^2) = 0
    v = roots(coeff_mat);
    
    v = v(imag(v)==0);       %% keeping only real volumes
    
    vmax = max(v);
    vmin = min(v); 
    func = @(v) [log(v(3)/(v(3)-b))-a*alpha*log((sigma*b+v(3))/(eps*b+v(3)))/((sigma-eps)*b*R*T)+log(1/v(3))+v(1)*v(3)/(R*T)-(log(v(2)/(v(2)-b))-a*alpha*log((sigma*b+v(2))/(eps*b+v(2)))/((sigma-eps)*b*R*T)+log(1/v(2))+v(1)*v(2)/(R*T));v(1)*v(2)^3+(v(1)*b-R*T)*v(2)^2+(a*alpha-2*R*T*b-3*v(1)*b^2)*v(2)+(R*T*b^2+v(1)*b^3-a*alpha*b);v(1)*v(3)^3+(v(1)*b-R*T)*v(3)^2+(a*alpha-2*R*T*b-3*v(1)*b^2)*v(3)+(R*T*b^2+v(1)*b^3-a*alpha*b)];
    xo = [Psat;vmax;vmin];
    [v,fval]= fsolve(func,xo);         %% v(1),v(2),v(3) respectively are pressure, vmax and vmin  
    P_array(i) = v(1);
    vg(i) = v(2);
    vl(i) = v(3);
    i = i+1;
end
x = vg(1):0.1:vg(20);
s = spline(vg,P_array,x);
plot(x,s,'b--o');
hold on;