U_i = 20;            % Ambient velocity
a = 4;               % cylinder radius

c = -a*5;            % starting coordinate (x)
b = a*10;            % ending coordinate (x)
d = -60;             % starting coordinate (y)
e = 60;              % ending coordinate (y)

n = a*50;            % number of intervals (step size in grid)

[x,y] = meshgrid([c:(b-c)/n:b],[d:(e-d)/n:e]');

for i = 1:length(x)
    for k = 1:length(x);
        f = sqrt(x(i,k).^2 + y(i,k).^2);
        if f < a
            x(i,k) = 0;
            y(i,k) = 0;
        end
    end
end

% Definition of polar variables
r = sqrt(x.^2+y.^2);
theta = atan2(y,x);                               

%% Creation of Streamline function 
z = U_i.*r.*(1-a^2./r.^2).*sin(theta);%- G*log(r)/(2*pi);

%% Creation of Figure

m = 100; 
s = ones(1,m+1)*a;
t = [0:2*pi/m:2*pi];

%% Streamline plot

contour(x,y,z,50);
hold on 
polar(t,s,'-k');
% axis square
title('Stream Lines');
grid off