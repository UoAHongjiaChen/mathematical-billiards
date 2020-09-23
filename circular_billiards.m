clear all; close all; clc; 

% Cobwebbing

% Our Difference Equation for Circular Billiards:
% x(n+1) = (x(n) + n*theta) mod2pi

% We test different values of theta to see if they exhibit periodic orbits

%% Period 4 orbit with lines and dots
circular_billiards_cobweb(0, pi/4, 100, 'o-'); 
% The line plot creates a nice closed triangle to indicate periodicty but
% also doesn't show the discontinuity. It may be hard to follow.


%% Period 4 orbit with dots
for i=[0:0.3:2*pi]   % a range of starting conditions
    circular_billiards_cobweb(i, pi/4, 100, 'o');
end

% The dot plot clearly shows the discotinuity and reflects the classical
% cobweb diagram. However we will need to manually cobweb to see the
% periodicity now.


%% Period 10 orbit with 3 complete rotations with lines
circular_billiards_cobweb(0, 3*pi/10, 100, 'o-'); 
% A nice closed 'graph' to indicate periodicity.


%% Period 10 orbit with 3 complete rotations with dots
for i=[0:0.3:2*pi]   % a range of starting conditions
    circular_billiards_cobweb(i, 3*pi/10, 100, 'o');  
end

% We need to manually cobweb this to find out periodicty.


%% Period 7 orbit with 4 complete roations with lines AND dots
circular_billiards_cobweb(0, 4*pi/7, 100, 'o-'); 
% A different closed 'graph' to indicate periodicity.


%% Non-periodic orbit with lines
circular_billiards_cobweb(0, pi/4 + 0.01, 1000, '-'); 
% Easy to see non-periodicty, if we were unsure we could increase the
% number of iterations. If the picture keeps changing, i.e: more and more
% area becomes filled, it is non-periodic.


%% Non-periodic orbit with dots
circular_billiards_cobweb(0, pi/4 + 0.01, 1000, 'o');
% Again we need to manually cobweb this to determine periodicity.


