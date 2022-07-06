clear;
clc;
syms w z x
u = 110;                             %free stream velocity
c = 1;                               %chord length
b = c/4;
ccmax = [0.03 0.05 0.06 0.07];       % max camber to chord ratios
tmaxc = [0.05 0.07 0.1 0.12];        %max thickness to chord ratios
alphaa = [2 4 9 11];                 %angles of attack in degrees
for i = 1:4
alpha = alphaa(i)*pi/180;
% circle parameters in w-plan
beta = 2*ccmax(i);                   
e = (tmaxc(i))/1.3;
a = b*(1+e)/cos(beta);
w0 = complex(-b*e,a*(beta));
gamma = 4*a*u*pi*sin(alpha+beta);    %circulation formula
gg = complex(0,gamma);
e_ialha = complex(cos(alpha),sin(alpha));  
e_n_ialpha = complex(cos(alpha),-sin(alpha));
F = (u*((w-w0)*e_n_ialpha+(e_ialha*a^2/(w-w0)))+(gg/(2*pi))*log((w-w0)*e_n_ialpha/a));   %complex potential function in w plan
z = w + b^2/w;                       % mapping w into z plan
dFdw = diff(F,w);
dzdw = diff(z,w);
dFdz = dFdw/dzdw;
cp = 1 - ((abs(dFdz))^2)/(u^2);      %pressure cooficient function in w
v_u = abs(dFdz)/u;                   %velocity magnitude function in w
jo = 0;
for theta = pi:-0.03:0               %upper surface from leading to trailing edge
jo = jo + 1;
r(jo) = b*(1+e*(1-cos(theta))+beta*sin(theta));   %r function in theta for w-plan circle
www(jo) = r(jo)*complex(cos(theta),sin(theta));   %value of w at certain theta
cpp(jo) = subs(cp,w,www(jo));                     %value of cpp at certain theta
v_u_all(jo) = subs(v_u,w,www(jo));                %magnitude of velocity at certain theta
airfoil(jo) = subs(z,w,www(jo));                  %point at the airfoil at certain theta
ree(jo) = real(airfoil(jo));                      %x-value of airfoil at certain theta
immag(jo) = imag(airfoil(jo));                    %y-value of airfoil at certain theta
end
figure(i)
l = -0.5:1/(length(cpp)-1):0.5;
plot(l,cpp,'LineWidth',1.6,'color','r')
title(sprintf('cp of airfoil for case %s',num2str(i)))
xlabel('x')
ylabel('cp')
hold on
plot(ree,immag,'LineWidth',1.6,'color','r')
hold on
figure(i+4)
plot(l,v_u_all,'LineWidth',1.6,'color','r')
title(sprintf('velocity to free stream velocity ratio of airfoil for case %s',num2str(i)))
xlabel('x')
ylabel('v/u')
hold on
cpp = [];
ree = [];
immag = [];
l = [];
jo = 0;
for theta = -pi:0.03:0               %lower surface from leading to trailing edge
jo = jo + 1;
r(jo) = b*(1+e*(1-cos(theta))+beta*sin(theta));   %r function in theta for w-plan circle
www(jo) = r(jo)*complex(cos(theta),sin(theta));   %value of w at certain theta
cpp(jo) = subs(cp,w,www(jo));                     %value of cpp at certain theta
v_u_all(jo) = subs(v_u,w,www(jo));                %magnitude of velocity at certain theta
airfoil(jo) = subs(z,w,www(jo));                  %point at the airfoil at certain theta
ree(jo) = real(airfoil(jo));                      %x-value of airfoil at certain theta
immag(jo) = imag(airfoil(jo));                    %y-value of airfoil at certain theta
end
figure(i)
l = -0.5:1/(length(cpp)-1):0.5;
plot(l,cpp,'LineWidth',1.6,'color','b')
hold on
plot(ree,immag,'LineWidth',1.6,'color','b')
legend('upper surface','','lower surface','','Location','southeast' )
hold off
figure(4+i)
plot(l,v_u_all,'LineWidth',1.6,'color','b')
hold off
legend('upper surface','lower surface','Location','northeast' )
cpp = [];
v_u_all = [];
ree = [];
immag = [];
l = [];
jo = 0;
end
