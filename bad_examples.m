function bad_examples()
% In the following examples, V(x) is NOT LLF 
    x=mpolyfun.singles(2);
    disp('-----Van Der Pol system------');
    f=sys_VanDerPol(); 
    vdp=agrasys(f);vdp_win=agragui(vdp);
    V=bad_fun();
    vdp.set_LF(V);
    vdp.disp()
    gamma=vdp.max_level()
    vdp_win.set_plane(1,2,[-5 5],[-5 5]);
    vdp_win.line_width(2);% line width
    vdp_win.plot_density(50,50);% discretization density for curve plotting
    vdp_win.see_plane();
    vdp_win.show_everything('rbg');
    
    disp('-----Ratschan system------');
    f=sys_ratschan();V=fun_ratschan();
    sys8=agrasys(f,V);sys8_win=agragui(sys8);
    sys8.disp()
    sys8.switch_LLF_check_OFF();
    gamma=sys8.max_level()
    sys8_win.set_plane(1,2,[-10 10],[-10 10]);
    sys8_win.line_width(2);
    sys8_win.plot_density(50,50);
    sys8_win.see_plane();
    sys8_win.show_everything('rbg');
    
    % Systems and corresponding Lyapunov functions
    
    function f=sys_VanDerPol()
    % Time reversed Van der Pol dynamics
        f(1)=-x(2);
        f(2)=x(1)+x(1)^2*x(2)-x(2);
    end

    function V=bad_fun()
    % One can obtain this function by calling lyap(A,eye(2))
             V=1.5*x(1)^2+x(1)*x(2)+x(2)^2;
    end
   % The following example taken from
   % Henning Burchardt, and Stefan Ratschan
   % Estimating the Region of Attraction of Ordinary Differential Equations by Quantified Constraint Solving
   % Proceedings of the 3rd WSEAS International Conference on DYNAMICAL SYSTEMS and CONTROL (CONTROL'07),
    % WSEAS Press, 2007, 241-246
    function f=sys_ratschan()
        f(1)=x(2)+0.2*x(1)-0.1*x(1)^3+0.01*x(1)^2*x(2)^2;
        f(2)=-x(1);
    end
     
    function V=fun_ratschan()
        V=7*x(1)^2-2*x(1)*x(2)+3*x(2)^2;
    end    

end