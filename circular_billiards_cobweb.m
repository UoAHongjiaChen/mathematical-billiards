% MATLAB needs the first function to be the same name as the filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We draw a plot on the x_{n} - x_{n+1} plane.

% Our Difference Equation for Circular Billiards:
% x(n+1) = (x(n) + 2*theta) mod2pi

% x0 = inital condition/position
% theta = angle (constant)
% N = number of iterations


% The types (symbols) of plots we can make have the iterates being:
% dots: 'o'
% connected lines: '-'
% dotted lines: 'o-'

% Use "o" to see the discontinuity and manual cobwebbing.
% Use "-" for easier visuals of periodicity.

function x = circular_billiards_cobweb(x0, theta, N, type)
    
x = zeros(1, N+1);
x(1) = x0;

for n=1:N
   x(n+1) =  mod(x(n) + 2*theta, 2*pi);
end

plot(x(1:(N-1)), x(2:(N)), type), hold on  % depends on the symbol type the user wants
%plot(x(1:(N-1)), repelem(2*pi,1,N-1), 'g--')
plot([0:0.01:2*pi], [0:0.01:2*pi], 'r--')
title('Cobwebbing Circular Billiards')
xlabel('x_n')
ylabel('x_{n+1}')
%hold off
set(gca,'XLim',[0 2*pi], 'YLim', [0 2*pi])
end



