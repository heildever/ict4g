function x = polar_to_rectangle(M, P)
i=sqrt(-1);
a=M*cos(P*pi/180);
b=M*i*sin(P*pi/180);
x=a+b;
end