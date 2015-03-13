classdef agragui < handle
% Visualization of LF and ROA
% Developer: Max Demenkov, Institute of Control Sciences, Moscow, May 2013 - March 2015 
% max.demenkov@gmail.com
%
% sys=agragui(obj) where obj is an agrasys object
% h=sys.window() opens simple GUI
% "Standalone" functions:
% sys.set_plane(i,j,[xi_min xi_max],[xj_min xj_max]) sets cross-section
%                      plane spanned by x_i and x_j variables
% sys.see_plane() Actually opens the plot window
% sys.show_everything('abc') Here a,b,c - color letters for plotting \dot V=0, V=gamma and minima
% sys.set_simulation(T,tol)  T - simulation time, tol- size of "zero box"
% sys.show_ROA()  Draws cross-section of ROA computed by numerical simulation
    properties
   
        main=[];
        x_i=1;y_i=2; % Plane coordinates
        plot_linewidth=3; % thickness of lines on plots
        marker_size=3;
        x_lim=[-5 5];y_lim=[-5 5]; % plot windows limits
        x_density=50;y_density=50; % grid density along x and y axis
        roa_x_density=20;roa_y_density=20;
        tol=1e-3; % practical zero
        xbox=[];% zero box
        T=50; % simulation time
    end
    
    methods
        % class constructor
        function  sys=agragui(obj)
        % Copy reference to an agrasys object
              if isa(obj,'agrasys')
                   sys.main=obj;
                   n=sys.main.n;
                   sys.xbox=sys.tol*[-ones(n,1),ones(n,1)];
              else
                   error('First argument must be agrasys object');
              end
        end
        
        function h=window(sys)
            global agra
            agra=sys;
            h=agrawin();
        end
        
        function set_simulation(sys,T,varargin)
        % T - simulation time, tol- size of zero box
            sys.T=T;n=sys.main.n;
            if nargin>1
                tol=varargin{1};
                sys.xbox=tol*[-ones(n,1),ones(n,1)];
            end
        end
        
        function set_plane(sys,i,j,varargin)
         % Set cross-section plane
             sys.x_i=i;sys.y_i=j;
             if nargin>3
                sys.x_lim=varargin{1};sys.y_lim=varargin{2};
             end
        end
        
        function see_plane(sys)
             h=figure;axis([sys.x_lim sys.y_lim]);box on;grid on;hold on;
        end
        
         function plot_density(sys,x_density,y_density)
             sys.x_density=x_density;sys.y_density=y_density;
         end
          
         function line_width(sys,width)
             sys.plot_linewidth=width;
         end
         
         function plot_points(sys,pts,col)
         % Put points on figure
              for j=1:size(pts,2)
                  plot(pts(sys.x_i,j),pts(sys.y_i,j),[col '*'],'MarkerSize',sys.marker_size);
              end
         end
         
         function show_V(sys,col)
         % Draw cross-section of Lyapunov function
             C=sys.main.gamma_max;
             if C>sys.tol && C<inf
                     level_curve(sys,sys.main.V,C,col);
             end
         end
         
         function show_dVdt(sys,col)
         % Draw cross-section of \dot V(x) =0
                     level_curve(sys,sys.main.dV,0,col);
         end
         
         function show_everything(sys,col)
         % Plot \dot V=0, V=gamma, and minima    
          show_dVdt(sys,col(1));show_V(sys,col(2));plot_critical(sys,col(3));
         end
         
         function plot_critical(sys,col)
                     lmin=sys.main.minima();
                     plot_points(sys,lmin,col);
         end

         function show_ROA(sys)
         % Draws cross-section of ROA computed by numerical simulation
             global agra_f;agra_f=sys.main.f;
             
             color='g';
             
             [map,xticks,yticks]=Build_ROA(sys);
             m=length(yticks);
             n=length(xticks);
             for i=2:m-1
                  for j=2:n-1
                      if map(i,j)
                        plot(xticks(j),yticks(i),[color 'o'],'LineWidth',2,...
                         'MarkerEdgeColor',color,...
                         'MarkerFaceColor',color,...
                         'MarkerSize',10);
                      end
                  end
            end
         end
         
         function show_trajectory(sys)
         % Pick initial conditions on the plane and plot a trajectory
         global agra_f; agra_f=sys.main.f;
         
               [x_1,x_2]=ginput(1);
               [x,y]=traj(sys,x_1,x_2);plot(x,y,'b');
         end
         
         function flag=in_ROA(sys,x0)
         % IF final state in the box THEN this point belongs to ROA
           WarningState=warning;warning off;
           [TOUT,YOUT] = ode45(@eval_sys,[0 sys.T],x0);
           warning(WarningState);
           x_end=YOUT(end,:)';
             if all(sys.xbox(:,1)<=x_end) && all(x_end<=sys.xbox(:,2))
                flag=1;
             else
                flag=0;
             end
         end
    
         function [x,y]=traj(sys,x0,y0)
         % Computes single trajectory
           WarningState=warning;warning off;
           n=sys.main.n;
           xvec=zeros(1,n);xvec(sys.x_i)=x0;xvec(sys.y_i)=y0;
           [TOUT,YOUT] = ode45(@eval_sys,[0 sys.T],xvec);
           warning(WarningState);
           x=YOUT(:,sys.x_i)';y=YOUT(:,sys.y_i)';
         end
    
         function [map,x,y]=Build_ROA(sys)
        % Build region of attraction in the form of 0/1 2D map
                     
                     n=length(sys.main.f);
                     [X,x,y]=xgrid(sys,sys.roa_x_density,sys.roa_y_density,n);
                     map=zeros(length(y),length(x));
                     % for the waitbar /begin
                     ProgressTitle='Computing...';
                     h = waitbar(0,ProgressTitle);drawnow;
                     binc=0.05;bcount=binc;
                     N=length(y)*length(x);
                     k=1;
                     % for the waitbar /end
                     for i=1:length(y)
                         for j=1:length(x)
                             x0=X(i,j,1:n);
                             map(i,j)=in_ROA(sys,x0);
                             % for the waitbar /begin
                              if (k/N)>=bcount
                                  waitbar(k/N);bcount=bcount+binc;
                              end
                              k=k+1;
                              % for the waitbar /end
                         end
                     end
                     close(h);
         end
       
        function level_curve(sys,fpoly,lval, col)
        % Draws level curves of fpoly(x)=lval
   
            n=varnum(fpoly);
            [X,x,y]=xgrid(sys,sys.x_density,sys.y_density,n);
            Z=fval(fpoly,X);
            contour(x,y,Z,col,'LevelList',lval,'LineWidth',sys.plot_linewidth);
     
       end % level_curve
     
      function  [X,x,y]=xgrid(sys,xd,yd,n)
      % Generate 2D grid of points; return 3D array 
                   dx=(sys.x_lim(2)-sys.x_lim(1))/xd;
                   dy=(sys.y_lim(2)-sys.y_lim(1))/yd;
                   x=sys.x_lim(1):dx:sys.x_lim(2);
                   y=sys.y_lim(1):dy:sys.y_lim(2);
                   xvec=zeros(1,n);
                   X=zeros(length(y),length(x),n);
                   for i=1:length(y)
                       for j=1:length(x)
                           xvec(sys.x_i)=x(j);xvec(sys.y_i)=y(i);
                           X(i,j,1:n)=xvec;
                       end
                  end
      end
      
    end
end

function dx=eval_sys(t,x)
global agra_f;
         dx=zeros(length(x),1);
         for i=1:length(x)
             dx(i)=fval(agra_f(i),x');
         end
end