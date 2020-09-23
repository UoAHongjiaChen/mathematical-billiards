% THIS HAS AN AREA FOR THE FORMULA OF A TRAPEZIUM FOR HITTING OPPOSITE SIDE
% YET IT WORKS PERFECTLY

% function y = square_area(alpha, P, period)
% % Computes the 'area' of a billiard trajectory within the square.
% % User specifies the initial angle 'alpha' and initial position 'P' as well
% % as the period of the periodic orbit that these initial conditions yield/
% 
% [alpha, P] = square_map(alpha, P, period + 1);  % plus one as it comes back around
% 
% side_jump = mod(diff(floor(P)),4);  % 1 for right adj, 2 for opp, 3 for left adj
% 
% area = 0;
% 
% for i=1:period
%     x_i = P(i) - floor(P(i)); 
%     x_i_p1 = P(i+1) - floor(P(i+1)); 
%     
%     % BELOW WORKS PERFECTLY FOR SOME REASON
%     if side_jump == 3
%         area = area + (1 - 0.5*x_i)*x_i_p1;
%     else
%         area = area + 0.5*(1-x_i)*x_i_p1;
%     end    
% 
% end
% y = period-area;  % could just have y = area and we our unique periodic orbit minimizes y
% end





function y = square_area(alpha, P, period)
% Computes the 'area' of a billiard trajectory within the square.
% User specifies the initial angle 'alpha' and initial position 'P' as well
% as the period of the periodic orbit that these initial conditions yield/

[alpha, P] = square_map(alpha, P, period + 1);  % plus one as it comes back around

side_jump = mod(diff(floor(P)),4);  % 1 for right adj, 2 for opp, 3 for left adj

area = 0;

for i=1:period
    x_i = P(i) - floor(P(i)); 
    x_i_p1 = P(i+1) - floor(P(i+1)); 

    if side_jump(i) == 3
        area = area + (1 - 0.5*x_i)*x_i_p1;
    elseif side_jump(i) == 2
        area = area + 0.5*(1-x_i+x_i_p1);
    elseif side_jump(i) == 1
        area = area + 0.5*(1-x_i)*x_i_p1;
    end    
end
y = period-area;  % could just have y = area and we our unique periodic orbit minimizes y
end

