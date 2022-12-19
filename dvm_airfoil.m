clc 
clear all
close all
nacaseries = input('Enter the 4-digit naca series = ','s');
n = input('Enter the number of nodes = ');
 c = input('Enter the chord length = ');
 U_inf = input('Enter the freestream velocity = ');
 rho = input('Enter the freestream density = ');
 alpha = input('Enter angle of attack (in degrees) = ');

 s1 = str2double(nacaseries(1));
 s2 = str2double(nacaseries(2));
 s3 = str2double(nacaseries(3));
 s4 = str2double(nacaseries(4));
    m = s1*0.01; p = s2*0.1 ; t = (10*s3+s4)*0.01;


for i= 1:n
    
    theta = (i-1)*2*pi/n;
    x = 0.5*c*(1+cos(theta));
if(x/c)<p
    yc(i) = m*c/p^2*(2*p*(x/c)-(x/c)^2);
    dydx(i) = (2*m/p^2)* (p-x/c);
    beta(i) = atan(dydx(i));
else
    yc(i) = m*c/(1-p)^2 * ((1-2*p)+2*p*(x/c)-(x/c)^2);
    dydx(i) = (2*m/(1-p)^2)* (p-x/c);
    beta(i) = atan(dydx(i));
end
yt=5*t*c*(0.2969*sqrt(x/c)-0.1260*(x/c)...
    -0.3516*(x/c)^2+0.2843*(x/c)^3-0.1036*(x/c)^4);

% plot(x,yc,'*r')
% hold on

if(i<(0.5*n+1))
    xa(i)=x - yt*sin(beta(i));
    ya(i)=yc(i)+yt*cos(beta(i));
else
    xa(i)=x + yt*sin(beta(i));
    ya(i)=yc(i)-yt*cos(beta(i));
end

end
xa(n+1)= c ; 
ya(n+1) = 0; 
yc(n+1) = 0;  % trailing edge
% plot(xa,ya,'k-')
hold on
%plot(xa,ya,'k -')

k = 0;
for i = 0.5*n+1:-1:1
    k = k+1;
    x_panel(k) = xa(i);
    y_panel(k) = ya(i);
end
for i = n:-1:0.5*n+1
    k = k+1;
    x_panel(k) = xa(i);
    y_panel(k) = ya(i);
end

 plot(x_panel,y_panel,'-*r')
axis equal
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)

% No. of panels and vortex placement

for i = 1:k-1
    dx(i) = x_panel(i+1)-x_panel(i);
    dy(i) = y_panel(i+1)-y_panel(i);
    s(i) = sqrt(dx(i)^2+dy(i)^2);     % length of panels
    nx(i) = -dy(i)/s(i);
    ny(i) = dx(i)/s(i);
end

% location of vortex and collocation point
for i = 1:k-1
    if i<=0.5*n+1
    x_vort(i) = x_panel(i) + 0.25*dx(i);
    y_vort(i) = y_panel(i) + 0.25*dy(i);
    x_coll(i) = x_panel(i) + 0.75*dx(i);
    y_coll(i) = y_panel(i) + 0.75*dy(i);
    else
    x_vort(i) = x_panel(i+1) - 0.25*dx(i);
    y_vort(i) = y_panel(i+1) - 0.25*dy(i);
    x_coll(i) = x_panel(i+1) - 0.75*dx(i);
    y_coll(i) = y_panel(i+1) - 0.75*dy(i);  
    end
end
plot(x_vort,y_vort,'og')
plot(x_coll,y_coll,'db')

G = zeros(k-1,k-1);

for i = 1:k-1
    for j = 1:k-1
    x_dist = x_coll(i)-x_vort(j);
    y_dist = y_coll(i)-y_vort(j);
    dist = x_dist^2+y_dist^2;
    num = y_dist*nx(i)-x_dist*ny(i);
    deno = 2*pi*dist;
    G(i,j) = num/deno;
    end
    R(i) = -U_inf*(cosd(alpha)*nx(i)+sind(alpha)*ny(i));
end

gamma = G^-1*R';

lift = 0;
for i = 1:k-1
    lift = lift+rho*U_inf*gamma(i);
end
drag = 0;

for i = 1:k-1
    for j = 1:k-1
        if i~=j
        x_dist = x_vort(i)-x_vort(j);
        y_dist = y_vort(i)-y_vort(j);
        dist = x_dist^2+y_dist^2;
        num = y_dist*sind(alpha)+x_dist*cosd(alpha);
        deno = 2*pi*dist;
        drag = drag+rho*gamma(i)*gamma(j)*num/deno;
        end
    end
end

dy_press = 0.5*rho*U_inf^2*c;
cl = lift/dy_press
cd = drag/dy_press
