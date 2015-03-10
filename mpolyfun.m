classdef mpolyfun
% Auxiliary functions for mpoly polynomials
% Developer: Max Demenkov, Institute of Control Sciences, Moscow, May 2013 - March 2015
% max.demenkov@gmail.com
%
% x=mpolyfun.singles(n) creates n-dim vector of monomials of one variable each
%                                  Example: x=mpolyfun.singles(2); f=x(1)*x(2);
% p=mpolyfun.dot_product(x1,x2) returns mpoly polynomial <x1,x2>
% G=mpolyfun.grad(f) returns either gradient (mpoly array) of polynomial function f  
%                                    or the Jacobian if an array f(1),f(2)...f(n) is provided
% H=mpolyfun.hessian(f)  returns Hessian (mpoly array) of polynomial function
% M=mpolyfun.mpoly2real(H,x) converts mpoly array H into real matrix M at point x
% A=mpolyfun.get_A(f)  extracts linearization matrix A from dynamical system f
% V=mpolyfun.get_quadratic(Q) translates positive definite matrix Q into
%                                        polynomial V=x'QX
% V=mpolyfun.get_lyapunov(f,C)  Solves Lyapunov equation A'Q+QA'=-C
%                                          where \dot x=Ax is the linearization of \dot x=f(x)
%                                         and returns polynomial V=x'Qx
%  Q=mpolyfun.get_approximation(V) Returns matrix Q of quadratic approximation
%                                                x'Qx of V(x) around x=0
%  P=mpolyfun.sdprojection(Q,tol)  Projection of matrix Q (in Frobenius norm) 
%                                                on the semidefinite cone
 
    methods(Static)
        
            function p=dot_product(g1,g2)
            % p=<g1,g2>
                p=mpoly;
                for i=1:length(g1)
                    p=p+g1(i)*g2(i);
                end
            end
            
            function G=grad(fun,varargin)
            % grad (f) returns gradient of polynomial function f
            % or Jacobian if length(f) >1
                 if nargin>1
                    n=varargin{1};
                 else
                    n=varnum(fun(1));
                 end
                 m=length(fun);
                 G(m,n)=mpoly;
                 for j=1:m
                     for i=1:n
                         G(j,i)=derivative(fun(j),i);
                     end
                 end
            end
         
            function H=hessian(fun,varargin)
            % returns Hessian of polynomial function
              if nargin>1
                    n=varargin{1};
                   else
                    n=varnum(fun);
              end
              H(n,n)=mpoly;g=mpolyfun.grad(fun,n);
              for i=1:n
                  for j=1:n
                      H(i,j)=derivative(g(j),i);
                  end
              end
            end
            
            function  x=singles(n)
            % create vector of monomials of one variable each
                         x(n)=mpoly;
                         for i=1:n
                             zvec=zeros(1,n);zvec(i)=1;
                             x(i)=mpoly(1,zvec);
                         end
            end
            
         function qpoly=get_quadratic(Q)
          % Translate positive definite matrix Q into polynomial
              n=size(Q,1);m=size(Q,2);
              if n~=m, error('Matrix is not square'); end
              if ~all(eig(Q)>=0), error('Matrix is not positive definite'); end
              p(n)=mpoly;qpoly=mpoly;
              x=mpolyfun.singles(n);
              for i=1:n
                  % matrix-vector multiplication p=Q*x
                  p(i)=mpolyfun.dot_product(Q(i,1:end),x);
              end
              % dot product x'*p
              qpoly=mpolyfun.dot_product(x,p);
         end
         
         function P=sdprojection(Q,tol)
         % Projection of matrix Q (in Frobenius norm) on the semidefinite cone
           [V,E]=eig(Q); % Obtain eigenvectors and eigenvalues
           d=diag(E);
           % If matrix is positive definite, do nothing
            if all(d>=tol), P=Q; return; end
            % Otherwise, set negative eigenvalues to tol>0
            d=max(d,tol);
            P=V*diag(d)*V';
         end
         
         function V=get_lyapunov(f,C)
         % Solves Lyapunov equation A'Q+QA'=-C
         % where \dot x=Ax is the linearization of \dot x=f(x)
         % and returns polynomial x'Qx
         % Matlab Control Toolbox must be installed
             A=mpolyfun.get_A(f);
             Q=lyap(A',C);
             V=mpolyfun.get_quadratic(Q);
         end
         
         function V=axes2ellipsoid(a,b,c)
         % Generates polynomial function from
         % the given ellipsoid semi-principal axes
         % f(x,y,z)=x^2/a^2+y^2/b^2+z^2/c^2
            x=mpolyfun.singles(3);
            V=(1/a^2)*x(1)^2+(1/b^2)*x(2)^2+(1/c^2)*x(3)^2;
         end
         
         function A=get_A(f)
         % Get linearization matrix A from dynamical system f
             n=varnum(f(1));
             J=mpolyfun.grad(f);% Get the Jacobian
             A=mpolyfun.mpoly2real(J,zeros(1,n));
         end
         
         function Q=get_approximation(V)
         % Return matrix Q of quadratic approximation
         % x'Qx of V(x) around x=0
              n=varnum(V);
              H=mpolyfun.hessian(V);
              Q=0.5*mpolyfun.mpoly2real(H,zeros(1,n));
         end
         
         function R=mpoly2real(M,x)
         % M is 2D array of mpoly objects
            [m,n]=size(M);
            R=zeros(m,n);
            for i=1:m
                for j=1:n
                    R(i,j)=fval(M(i,j),x);
                end
            end
         end
         
         function flag=is_LLF(f,V)
         % Check if LLF is suitable for the given system
             A=mpolyfun.get_A(f);% Get linearization of the system
             Q=mpolyfun.get_approximation(V); % Get quadratic approximation of the LF
             P=(A')*Q+Q*A;% Check Lyapunov inequality
             if all(eig(P)<0)
               flag=true;
             else
               flag=false;
             end
         end
         
    end
end    