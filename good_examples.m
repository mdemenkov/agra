function good_examples()
% Systems and corresponding Lyapunov functions
% taken from S. Ratschan, Z. She, Providing a basin of attraction to a target region 
% of polynomial systems by computation of Lyapunov-like functions.
% SIAM Journal of Control and Optimization, Vol. 48, No. 7,  pp. 4377-4394, 2010.

    x=mpolyfun.singles(2);
    
    disp('-----Van Der Pol system------');
      % Time reversed Van der Pol dynamics
        f(1)=-x(2);
        f(2)=x(1)+x(1)^2*x(2)-x(2);
        a=1;b=-0.3449172;c=0.8589766;
        V=mpoly([a;b;c],[2 0;1 1;0 2]);
        
    vdp=agrasys(f,V);vdp_win=agragui(vdp);
    vdp.disp()
    gamma=vdp.max_level()

    minima=vdp.minima()
    vdp_win.set_plane(1,2,[-2 2],[-2 2]);
    vdp_win.line_width(3);% line width
    vdp_win.see_plane();
    vdp_win.show_everything('rbg');
    V=mpolyfun.get_lyapunov(f,eye(2)); % LF from Lyapunov equation
    vdp.set_LF(V);
    vdp.max_level();
    vdp_win.line_width(1);
    vdp_win.show_everything('rbg');

    disp('-----SOS example------------');
       f(1)=-x(1)+x(2);
       f(2)=0.1*x(1)-2*x(2)-x(1)^2-0.1*x(1)^3;
       % LF is always positive
       V=x(1)^2+1.6513*x(2)^2;
       
    sos=agrasys(f,V);sos_win=agragui(sos);
    sos
    gamma=sos.max_level()
    minima=sos.minima()
     sos_win.set_plane(1,2,[-10 4],[-10 4]);
     sos_win.plot_density(50,50);% discretization density for curve plotting 
     sos_win.line_width(3);
     sos_win.see_plane();
     sos_win.show_everything('rbg');

     % Now, lets change LF and see what happens
     % A little bit artificial example
     P=mpolyfun.get_approximation(V);
     sos_win.line_width(1);
     for i=1:10
         P(1,1)=P(1,1)+0.5;% change first coefficient of LF
         Q=mpolyfun.sdprojection(P,1e-3); % project coefficients of (candidate) LF on the "semidefinite set" of coefficients
         V=mpolyfun.get_quadratic(Q);% convert to mpoly
         sos.set_LF(V);sos.max_level();
         sos_win.show_dVdt('r');sos_win.show_V('b');
     end
 
end
