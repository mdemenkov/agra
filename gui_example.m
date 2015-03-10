function gui_example()
% Simple GUI demonstration
       x=mpolyfun.singles(2);
       f(1)=-x(1)+x(2);
       f(2)=0.1*x(1)-2*x(2)-x(1)^2-0.1*x(1)^3;
       V=x(1)^2+1.6513*x(2)^2;
       sos=agrasys(f,V);
       sos.max_level();

       sos_gui=agragui(sos);
       sos_gui.set_simulation(30,0.05);
       sos_gui.window();
end