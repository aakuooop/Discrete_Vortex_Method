clc 
clear all
close all
nacaseries = input('Enter the 4-digit naca series = ','s');
 c = input('Enter the chord length = ');
  U_inf = input('Enter the freestream velocity = ');
 rho = input('Enter the freestream density = ');
 alpha = input('Enter the angle of attack (in degrees) = ');
 h_qcp = input('Enter the height of quarter-chord point from ground (Ground clearance) = ');
 
 h = h_qcp-(0.75*c*sind(alpha)); 
 

 s1 = str2double(nacaseries(1));
 s2 = str2double(nacaseries(2));
 s3 = str2double(nacaseries(3));
 s4 = str2double(nacaseries(4));
    m = s1*0.01; p = s2*0.1 ; t = (10*s3+s4)*0.01;

n = 50;
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

% at height 'h' inclined at AOA = alpha


for i = 1:n+1
   xa(i) = xa(i)-c;
end
 
alpha = -alpha;

% rotating airfoil by angle 'alpha' about TE
for i = 1:n+1
    xa(i) = xa(i)*cosd(alpha)-ya(i)*sind(alpha);
    ya(i) = xa(i)*sind(alpha)+ya(i)*cosd(alpha);
    yc(i) = xa(i)*sind(alpha)+yc(i)*cosd(alpha);
end

for i = 1:n+1
    ya(i) = ya(i)+h;
    yc(i) = yc(i)+h;
    xa(i) = xa(i)+c;
end

k = 0;
for i = 0.5*n+1:-1:1
    k = k+1;
    x_panel(k) = xa(i);
    y_panel(k) = yc(i);
    x_imag(k) = xa(i);
    y_imag(k) = -yc(i);
end
plot(xa,ya,'-k')
hold on
 plot(x_panel,y_panel,'-*r')
 plot(xa,-ya,'-k')
  plot(x_imag,y_imag,'-*r')
 axis equal
 xlabel('x','fontsize',15)
 ylabel('y','fontsize',15)

% No. of panels and vortex placement

for i = 1:n/2
    dx(i) = x_panel(i+1)-x_panel(i);
    dy(i) = y_panel(i+1)-y_panel(i);
    s(i) = sqrt(dx(i)^2+dy(i)^2);     % length of panels
    nx(i) = -dy(i)/s(i);
    ny(i) = dx(i)/s(i);
end

% location of vortex and collocation point 
for i = 1:n/2
    x_vort(i) = x_panel(i) + 0.25*dx(i);
    y_vort(i) = y_panel(i) + 0.25*dy(i);
    x_coll(i) = x_panel(i) + 0.75*dx(i);
    y_coll(i) = y_panel(i) + 0.75*dy(i);
end

% location of vortex point of image 

for i = 1:n/2
    x_vort_img(i) = x_vort(i);
    y_vort_img(i) = -y_vort(i);
end

plot(x_vort,y_vort,'og')
plot(x_coll,y_coll,'db')
plot(x_vort_img,y_vort_img,'og')
v1 = [-1,2];
v2 = [0,0];
line(v1,v2,'LineWidth',2)

G = zeros(n/2,n/2);

for i = 1:n/2
    for j = 1:n/2
    x_dist1 = x_coll(i)-x_vort(j);
    y_dist1 = y_coll(i)-y_vort(j);
    dist1 = x_dist1^2+y_dist1^2;
    num1 = y_dist1*nx(i)-x_dist1*ny(i);
    deno1 = 2*pi*dist1;
    term_1 = num1/deno1;
    
    % due to image on ith panel
    
    x_dist2 = x_coll(i)-x_vort_img(j);
    y_dist2 = y_coll(i)-y_vort_img(j);
    dist2 = x_dist2^2+y_dist2^2;
    num2 = y_dist2*nx(i)-x_dist2*ny(i);
    deno2 = -2*pi*dist2;
    term_2 = num2/deno2;
    
    G(i,j) = term_1+term_2;
    end
    R(i) = -U_inf*nx(i);
end

gamma = G^-1*R';

lift = 0;
drag = 0;
for i = 1:n/2
    lift_panel = 0;
    drag_panel = 0;
    lift_img = 0;
    drag_img = 0;
    for j = 1:n/2
        if i~=j
    x_dist1 = x_vort(i)-x_vort(j);
    y_dist1 = y_vort(i)-y_vort(j);
    dist1 = x_dist1^2+y_dist1^2;
    deno1 = 2*pi*dist1;
    lift_term_1 = (rho*gamma(i)*gamma(j)*y_dist1)/deno1;
    lift_panel = lift_panel+lift_term_1;
    drag_term_1 = (gamma(j)*x_dist1)/deno1;
    drag_panel = drag_panel+drag_term_1;
      end
    % due to image on ith panel
    
    x_dist2 = x_vort(i)-x_vort_img(j);
    y_dist2 = y_vort(i)-y_vort_img(j);
    dist2 = x_dist2^2+y_dist2^2;
    deno2 = 2*pi*dist2;
    lift_term_2 = (rho*gamma(i)*-gamma(j)*y_dist2)/deno2;
    lift_img = lift_img+lift_term_2;
    drag_term_2 = (-gamma(j)*x_dist2)/deno2;
    drag_img = drag_img+drag_term_2;
    
    lift_tot = lift_panel+lift_img;
    drag_tot = drag_panel+drag_img;
    
    end
    % lift and drag 
    lift_self = rho*U_inf*gamma(i);
    lift = lift+lift_self+lift_tot;
    drag = drag+rho*gamma(i)*drag_tot;
end

dy_press = 0.5*rho*U_inf^2*c;
cl = lift/dy_press
cd = drag/dy_press
