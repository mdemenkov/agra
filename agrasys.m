classdef agrasys < handle
% ROA estimation via numerical algebraic geometry
% PHCpack (http://homepages.math.uic.edu/~jan/download.html) must be installed
% Developer: Max Demenkov, Institute of Control Sciences, Moscow, May 2013 - March 2015 
% max.demenkov@gmail.com
%
%  sys=agrasys(f,V)  where f is dynamical system and V is local Lyapunov function
%  Here is a simple example:
%      x=mpolyfun.singles(2);
%      f(1)=-x(1)+x(2);
%      f(2)=0.1*x(1)-2*x(2)-x(1)^2-0.1*x(1)^3;
%      V=x(1)^2+1.6513*x(2)^2;
%
%  C=sys.max_level() finds maximum level of LLF
%  sys.minima() returns array of minimum points
% 
    properties
        f=[];% dynamical system
        n=0;% its dimension
        V=[];% Lyapunov function
        dV=[]; % Its derivative
        grad_V=[]; % gradient of LF
        grad_dV=[]; % gradient of \dot V
        gamma_max=0; % ROA estimation is given by { x | V(x)<=gamma_max }
        l_ext=[]; % Local extremums of Lyapunov function on dV/dt=0
        gmin=[]; % Global minima (might be several points)
        eqpoints=[]; % Equilibrium points
        tol=1e-6; % "practical zero"
        trueLLF=false;
        check_LLF=true;
        KKT=[];% KKT system
    end
    
    methods
        
        % class constructor
        function sys=agrasys(varargin)
            if nargin>0
               if length(varargin{1})<2
                  error('System dimension<2');
               end
               sys.f=varargin{1};
               sys.n=varnum(sys.f(1));
               if sys.n~=length(sys.f)
                  error('Number of equations is not equal state space dimension');
               end
               PHC_setup();
               if nargin>1
                   V=varargin{2};
                   sys.set_LF(V);
               end
            end
        end
        
        function set_LF(sys,V)
             % Changes candidate LF for the system
             grad_V=mpolyfun.grad(V);
             dV=mpolyfun.dot_product(sys.f,grad_V);
             grad_dV=mpolyfun.grad(dV);
             sys.V=V;
             sys.dV=dV;
             sys.grad_V=grad_V;
             sys.grad_dV=grad_dV;
             % Lets check if its suitable
             sys.trueLLF=mpolyfun.is_LLF(sys.f,V);
             sys.KKT=KKT_system(sys);
        end
        
        function switch_LLF_check_OFF(sys)
         % Switch off Lyapunov equation check
            sys.check_LLF=false;
            disp('LLF Hessian check is OFF');
        end
        
        function C=max_level(sys)
        % Finds maximum level of LLF
                   if ~sys.trueLLF
                       C=0;sys.gamma_max=0;
                       disp('LLF is bad, ROA=0');
                       if sys.check_LLF, return; end
                   else
                       disp('LLF is good');
                   end
                   sol=NAG_GetReal(sys.KKT,sys.tol);
                   if isempty(sol)
                     disp('No real solutions - path jumping or global stability');
                     C=inf;sys.gamma_max=inf;
                     sys.gmin=[];
                  else
                    % Separate x and Lagrange multipliers
                    rsol=sol(:,1:end-1);lambda=sol(:,end);
                    lvalues=fval(sys.V,rsol);
                    C=min(lvalues);
                    sys.gamma_max=C;sys.l_ext=rsol';
                    % Select all minimum points
                    ind=abs(lvalues-C)<=sys.tol;
                    nmin=sum(ind);
                    disp(['N of minimums=' num2str(nmin)]);
                    gmin=rsol(find(ind),:);lmin=lambda(find(ind))
                    Leval=zeros(nmin,1);
                    for i=1:nmin
                        % Check if in all minimum points gradients of V(x) and
                        % dV(x)/dt are pointing in the same direction
                        if lmin(i)>0, Leval(i)=1; end
                    end
                    if sum(Leval)<nmin
                        disp('Local minimum condition is not valid');
                        sys.gmin=[];
                    else
                        sys.gmin=gmin';
                    end
                  end
        end
         
          function eq=KKT_system(sys)
          % Return KKT system of polynomial equations:
          % grad(V)-\lambda*grad(\dot V)=0
          % \dot V=0
             nvars=varnum(sys.V);
             eq=lagrange(sys.grad_V,sys.grad_dV);
             eq(nvars+1)=extend(sys.dV,1);
          end
          
        function  eqpoints=equilibria(sys)
                      eqpoints=NAG_GetReal(sys.f,sys.tol)';
                      sys.eqpoints=eqpoints;
        end
        
        function mpoints=minima(sys)
                    mpoints=sys.gmin;
        end
        
        function disp(sys)
                     sys.disp_sys();
                     disp(' ');
                     sys.disp_lyap();
        end
        
         function disp_sys(sys)
           for i=1:length(sys.f)
               disp(['dx' num2str(i) '/dt=' p2str(sys.f(i))]);
           end
         end
        
         function disp_lyap(sys)
             if ~isempty(sys.V)
                disp(['V(x)=' p2str(sys.V)]);
             end
         end
         
end
    
         
end

% Auxiliary functions
         
           function l=lagrange(V,dV)
          % Return the (array of) Lagrange function
                nvars=varnum(V(1));
                m=length(V);
                lambda=mpoly(1,[zeros(1,nvars) 1]);
                l(m)=mpoly;
                for i=1:m
                     l(i)=extend(V(i),1)-lambda*extend(dV(i),1);
                end
           end
         
          function RSol=NAG_GetReal(eq,mytol)
          % Convert polynomial system to PHCpack format and call it to solve the system
          % Return real solutions

             MatForm=[];RSol=[];
             for i=1:length(eq)
                 cvec=eq(i).cvec;pvec=eq(i).pvec;
                 n=varnum(eq(i));
                 MatForm=[MatForm;[cvec pvec;zeros(1,n+1)]];
             end
             disp('Entering PHC...');drawnow update;
             s=solve_system(MatForm);
             ns=size(s,2);
             for i=1:ns
                 x=PHC2Mat(s(i),n);
                 if all(abs(imag(x))<=mytol)
                     RSol=[RSol;real(x)'];
                 end
             end
             disp(['Solved ! All=' num2str(ns) ', real=' num2str(size(RSol,1))]);
          end
         
          function x=PHC2Mat(sol,n)
          % Convert solution from PHC format 
             x=zeros(n,1);
             for i=1:n
                 x(i)=getfield(sol,['x' num2str(i)]);
             end
          end
    
          function PHC_setup()
             % Check PHClab
             win_flag=~isempty(findstr(computer,'PCWIN'));
             if exist('solve_system.m')==0
                phclab_name=find_here('PHClab',1);
                if ~isempty(phclab_name)
                   if win_flag
                      plab=[pwd '\' phclab_name];
                   else
                      plab=[pwd '/' phclab_name];
                   end
                      path(path,plab);
                else
                   error('Matlab interface PHClab not found');
                end
             end
            % Check the executable
            if win_flag, pname='phc.exe';else, pname='phc'; end
            if exist(pname)==2 
               set_phcpath([pwd '/' pname]);
            else
               error('PHC not found');
            end
           
          end

          function name=find_here(str,isdir)
                      curdir=dir;
                      name=[];
                      for i=1:length(curdir)
                          if ~isempty(findstr(curdir(i).name,str))
                              if curdir(i).isdir==isdir
                                 name=curdir(i).name;
                                 return;
                              end
                          end
                      end
          end


