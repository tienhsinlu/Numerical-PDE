%This code is based on the code from Nick Trefethen
clc;
%RK, AF, AM, BDF
method='AF';


switch method
    case 'RK'
        
          h=zeros(6,1);
          h(1)=plot([-8 8],[0 0]);
          hold on
          h(2)=plot([0 0],[-8 8]);
          w = 0; W = w; z = exp(1i*pi*(0:200)/100);
          for i = 2:length(z)                  % order 1
            w = w-(1+w-z(i)); 
            W = [W; w]; 
          end
          h(3)=plot(W, 'LineWidth', 2, 'Displayname', 'order 1');
          
          w = 0; W = w; 
          for i = 2:length(z)                  % order 2
            w = w-(1+w+.5*w^2-z(i)^2)/(1+w); 
            W = [W; w];
          end
          h(4)=plot(W, 'LineWidth', 2, 'Displayname', 'order 2');
          
          w = 0; W = w; 
          for i = 2:length(z)                  % order 3
            w = w-(1+w+.5*w^2+w^3/6-z(i)^3)/(1+w+w^2/2); 
            W = [W; w];
          end
          h(5)=plot(W, 'LineWidth', 2, 'Displayname', 'order 3');
          w = 0; W = w; 
          for i = 2:length(z)                  % order 4
            w = w-(1+w+.5*w^2+w^3/6+w.^4/24-z(i)^4)/(1+w+w^2/2+w.^3/6);
            W = [W; w]; 
          end
          h(6)=plot(W, 'LineWidth', 2, 'Displayname', 'order 4');
          axis([-5 2 -3.5 3.5]), axis square, grid on, 
          title('Stability Regions for Runge-Kutta Methods')
          legend(h(3:end), 'Location', 'best')
          
          
    case 'AF'
          h=zeros(5,1);
          h(1)=plot([-8 8],[0 0]);
          hold on
          h(2)=plot([0 0],[-8 8]);
          z = exp(1i*pi*(0:200)/100); r = z-1;
          s = 1; 
          h(3)=plot(r./s, 'b-', 'LineWidth', 2, 'Displayname', 'order 1');                  % order 1
          s = (3-1./z)/2; 
          h(4)=plot(r./s, 'r-', 'LineWidth', 2, 'Displayname', 'order 2');                  % order 2
          s = (23-16./z+5./z.^2)/12; 
          h(5)=plot(r./s, 'k-', 'LineWidth', 2, 'Displayname', 'order 3');                  % order 3
          axis([-2.5 .5 -1.5 1.5]), axis square, grid on
          title Adams-Bashforth
          legend(h(3:end), 'Location', 'best')   
          
          
    
    case 'AM'
          z = exp(1i*pi*(0:200)/100); r = z-1;     
          h=zeros(6, 1);
          h(1)=plot([-8 8],[0 0]);
          hold on
          h(2)=plot([0 0],[-8 8]);
          s = (5*z+8-1./z)/12; 
          h(3)=plot(r./s, 'b-', 'LineWidth', 2, 'Displayname', 'order 3');                    % order 3
          s = (9*z+19-5./z+1./z.^2)/24; 
          h(4)=plot(r./s, 'r-', 'LineWidth', 2, 'Displayname', 'order 4');           % order 4 
          s = (251*z+646-264./z+106./z.^2-19./z.^3)/720; 
          h(5)=plot(r./s, 'k-', 'LineWidth', 2, 'Displayname', 'order 5');    % 5
          d = 1-1./z;
          s = 1-d/2-d.^2/12-d.^3/24-19*d.^4/720-3*d.^5/160; 
          h(6)=plot(d./s, 'LineWidth', 2, 'Displayname', 'order 6'); % 6
          axis([-7 1 -4 4]), axis square, grid on, title Adams-Moulton
          legend(h(3:end), 'Location', 'best')
         
    case 'BDF'
          h=zeros(7, 1);
          h(1)=plot([-10 20],[0 0]);
          hold on
          h(2)=plot([0 0],[-15 15]);
          r = 0; 
          for i = 1:5
              r = r+(d.^i)/i; 
              h(i+2)=plot(r, 'LineWidth', 2, 'Displayname', "order "+num2str(i));
              hold on
          end   % orders 1-5
          axis([-10 20 -15 15])
          axis square
          grid on
          title('Region of Stability of BDF Methods (exterior of the curves)')
          legend(h(3:end), 'Location', 'best')
          
          
end