classdef mpoly
   % Implement multivariate polynomial class
   % Developer: Max Demenkov, Institute of Control Sciences, Moscow, May 2013 - July 2015
   % max.demenkov@gmail.com
   %

    properties
        cvec % Vector of monomial coefficients
        pvec % Matrix where each row contains powers for monomial variables
        % Example: If p(x)=c_1*x_1*x_2+c_2*x_2^2 then
        %          cvec=[c_1;c_2];pvec=[1 1;0 2];
    end
    
    methods(Static=true)
        
        function [cvec,pvec]=dblentries(cvec1,pvec1,cvec2,pvec2)
        % basic procedure for adding polynomials
        cvec=cvec1;pvec=pvec1;
        cvecr=cvec2;pvecr=pvec2;
        for i=1:length(cvec)
             remlist=[];
             for j=1:length(cvecr)
                 if all(pvec(i,:)==pvecr(j,:))
                    cvec(i)=cvec(i)+cvecr(j);
                 else
                    remlist=[remlist;j];
                 end
             end
             pvecr=pvecr(remlist,:);cvecr=cvecr(remlist);
         end
         if ~isempty(cvecr)
             cvec=[cvec;cvecr];pvec=[pvec;pvecr];
         end
        end
        
        function  [cvec_out,pvec_out]=rmvsimilar(cvec_in,pvec_in)
         % removes similar monomials
         cvec_out=[];pvec_out=[];
         cvec=cvec_in;pvec=pvec_in;
          for i=1:length(cvec_in)
              r=pvec_in(i,:);
              if isempty(mpoly.findpowers(pvec_out,r))
                 ind=mpoly.findpowers(pvec_in,r);
                 pvec_out=[pvec_out;r];
                 cvec_out=[cvec_out;sum(cvec_in(ind))];
              end
          end
        end
        
         function ind=findpowers(P,r)
                        [n,m]=size(P);
                        ind=[];
                        for k=1:n
                            if all(P(k,:)==r)
                               ind=[ind;k];
                            end
                        end
        end
            
        function  [cvec,pvec]=rmvzeros(cvec1,pvec1)
                  rmconst=1e-10;% removes monomials less than that
                  remlist=[];
                  for i=1:length(cvec1)
                      if abs(cvec1(i))>rmconst
                         remlist=[remlist;i];
                      end
                  end
                  cvec=cvec1(remlist);pvec=pvec1(remlist,:);
        end
      
    end
    
    methods
        function obj=mpoly(varargin)
                 switch nargin
                     case 0
                        obj.cvec=[];obj.pvec=[];
                     case 2
                         if any(varargin{2}<0)
                            error('Negative powers');
                         elseif ~all(round(varargin{2})==varargin{2})
                            error('Powers should be integers');
                         else
                         [obj.cvec,obj.pvec]=mpoly.rmvzeros(varargin{1},...
                                                            varargin{2});
                         end
                     otherwise
                         error('Number of inputs');
                 end
        end
        
        function display(obj)
                 disp(p2str(obj));
        end
        
        function wholestr=p2str(obj)
            cvec=obj.cvec;pvec=obj.pvec;
            wholestr=[];
            for i=1:length(cvec)
                if abs(cvec(i))==1 && any(pvec(i,:))
                    if sign(cvec(i))==-1 && i==1
                       str='-';
                    else
                       str='';
                    end
                else
                   if i>1
                      str=num2str(abs(cvec(i)));
                   else 
                      str=num2str(cvec(i));
                   end
                   if any(pvec(i,:))
                      str=[str '*'];
                   end
                end
                for j=1:size(pvec,2)
                    if pvec(i,j)>0
                       str=[str 'x' num2str(j)];
                       if pvec(i,j)>1
                          str=[str '^' num2str(pvec(i,j))];
                       end
                       if j<size(pvec,2) && any(pvec(i,j+1:end))>0 
                          str=[str '*'];
                       end
                    end
                end
                if i<length(cvec) 
                   if sign(cvec(i+1))<0, sgn='-'; else sgn='+'; end
                   str=[str sgn];
                end
                wholestr=[wholestr, str];
            end
        end
        
        function r = plus(obj1,obj2)
         % PLUS  Implement obj1 + obj2 
        if isa(obj1,'double')
             num.cvec=obj1;num.pvec=zeros(1,size(obj2.pvec,2));
             [cvec,pvec]=mpoly.dblentries(num.cvec,num.pvec,...
                                          obj2.cvec,obj2.pvec);
        elseif isa(obj2,'double')
             num.cvec=obj2;num.pvec=zeros(1,size(obj1.pvec,2));
             [cvec,pvec]=mpoly.dblentries(num.cvec,num.pvec,...
                                          obj1.cvec,obj1.pvec);
        else
         [cvec,pvec]=mpoly.dblentries(obj1.cvec,obj1.pvec,...
                                      obj2.cvec,obj2.pvec);
        end
        r=mpoly(cvec,pvec);
        end % plus
      
        function r=uminus(obj)
            r=mpoly(-obj.cvec,obj.pvec);
        end
        
      function r = minus(obj1,obj2)
         % MINUS Implement obj1 - obj2
         r=plus(obj1,uminus(obj2));
      end % minus
      
      function r=mpower(obj,p)
      % Raises polynomial into power p 
          if isa(obj,'mpoly') && isa(p,'double') && fix(p)==p && p>0
              % power must be positive integer
             if length(obj.cvec)>1
                 r=obj;
                 for i=1:p-1
                     r=mtimes(r,obj);
                 end
             else % single monomial
                 r=mpoly(obj.cvec,obj.pvec.*p);
             end
          else
                 error('Power error');
          end
      end
      
      function r = mtimes(obj1,obj2)
        % MTIMES   Implement obj1 * obj2 
        if isa(obj1,'double')
             cvec=obj2.cvec;pvec=obj2.pvec;
             r=mpoly(cvec*obj1,pvec);
        elseif isa(obj2,'double')
             cvec=obj1.cvec;pvec=obj1.pvec;
             r=mpoly(cvec*obj2,pvec);
        else
            cvec1=obj1.cvec;pvec1=obj1.pvec;
            cvec2=obj2.cvec;pvec2=obj2.pvec;
            n=size(pvec1,2);m=length(cvec1)*length(cvec2);
            cvec3=zeros(m,1);pvec3=zeros(m,n);
            k=0;
            for i=1:length(cvec1)
                for j=1:length(cvec2)
                    k=k+1;
                    cvec3(k)=cvec1(i)*cvec2(j);
                    pvec3(k,:)=pvec1(i,:)+pvec2(j,:);
                end
            end
            [cvec_out,pvec_out]=mpoly.rmvsimilar(cvec3,pvec3);
            r=mpoly(cvec_out,pvec_out);
        end
      end % mtimes
      
      function r=derivative(obj,nx)
          cvec=obj.cvec;pvec=obj.pvec;
          remlist=[];
          for i=1:length(cvec)
              if pvec(i,nx)>0
                 cvec(i)=cvec(i)*pvec(i,nx);
                 pvec(i,nx)=pvec(i,nx)-1;
                 remlist=[remlist;i];
              end
          end
          if ~isempty(remlist)
              r=mpoly(cvec(remlist),pvec(remlist,:));
          else
              r=mpoly();
          end
      end
      
      function r=extend(obj,n)
          cvec=obj.cvec;pvec=obj.pvec;
          pvec=[pvec zeros(length(cvec),n)];
          r=mpoly(cvec,pvec);
      end
      
      function Z=fval(obj,X)
               % return values of polynomial at points X
               cvec=obj.cvec;pvec=obj.pvec;
               nx=size(pvec,2);
               w=zeros(length(cvec),1);
               if ismatrix(X)
                  % matrix with points as its rows 
                  m=size(X,1); Z=zeros(m,1); 
                  for i=1:m
                       for j=1:length(cvec)
                           w(j)=cvec(j)*prod(X(i,:).^pvec(j,:));
                       end
                       Z(i)=sum(w);
                  end
               else
                  % 3D array, each 2D entry is a point
                  m=size(X,1);n=size(X,2);
                  Z=zeros(m,n);
                  for i=1:m
                      for j=1:n
                          x=reshape(X(i,j,1:nx),[1 nx]);
                          for k=1:length(cvec)
                              w(k)=cvec(k)*prod(x.^pvec(k,:));
                          end
                          Z(i,j)=sum(w);
                      end
                  end
               end
      end
      
      function n=varnum(obj)
               n=size(obj.pvec,2);
      end 
      
    end
    
end
