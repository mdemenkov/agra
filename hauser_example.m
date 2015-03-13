function hauser_example()
% From: J. Hauser, M. C. Lai, Estimating quadratic stability domains by nonsmooth optimization. 
% Proc. of the ACC, pp. 571-576,  Chicago, 1992

x=mpolyfun.singles(3);
f(1)=x(2)+2*x(2)*x(3);
f(2)=x(3);
f(3)=-0.5*x(1)-2*x(2)-x(3);

V=mpolyfun.get_lyapunov(f,eye(3)); % LF from Lyapunov equation
hauser=agrasys(f,V);
hauser.max_level()

mygui=agragui(hauser);
mygui.set_simulation(30,0.05);
mygui.set_plane(2,3,[-1 1],[-1 1]);
mygui.see_plane();xlabel('x_2');ylabel('x_3');
mygui.show_ROA();
mygui.show_everything('rbm');

end