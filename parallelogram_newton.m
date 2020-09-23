%% Area testing

h = 1; gamma = pi/2;

x = linspace(0.01, 1-0.01, 1000);   % Initial positions on base side

period = 10;
alpha = atan(4);     %atan(h/(1+h/tan(gamma)));

for j=1:length(x)
    y(j) = parallelogram_area(h, gamma, alpha, x(j), period);
 
end

plot(x,y)
xlabel('P_0'); ylabel('Area')


%% Maximum area testing
clear all; clc

h = 3.5; gamma = pi/2-0.01;

alpha0 = atan(h*3/2);
period = 10;


x = linspace(0.01, 1-0.01, 100);   % Initial positions on base side

for j=1:length(x)
    y(j) = parallelogram_area(h, gamma, alpha0, x(j), period);
end

plot(x,y)
xlabel('P_0'); ylabel('Area')

%%

%% Parallelogram Newton solver

clear all; clc
eps = 1e-6;   % for finite differencing

% Size of parallelogram
h = 3.16; gamma = pi/2;

period = 10; % what period are we seaching for?

alpha = 0.98; P = 0.49;   % Initial guess

guesses = [alpha P]; % First column: init_angle and init_pos  

for j=1:7
    guess = guesses(j,:)';  % column vector
    
    % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
    [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
    dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

    guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

    alpha = guess(1); P = guess(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP SOLUTION for when we hop out of range
    alpha = mod(alpha, pi); P = mod(P, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    guesses(j+1,:) = [alpha P];
    
    if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6   % relative error
        break
    end
end

guesses




%% Following a period-10 orbit by keeping gamma constant and varying the height.

clear all; clc
tic
period = 10; % what period are we seaching for?

gamma = pi/2; % keep this fixed for now

% We know for square that the angle and position that is a period-4 orbit
% and that maximises the area is angle=pi/4, pos = 0.5. This is our
% initial guess.

eps = 1e-6;   % for finite differencing

height_perturb = 1e-2; % perturbing height of rectangle

heights = [1];             % The heights, filled in for square already
alpha_P_mat = [atan(4/1) 0.125]; % Already filled in for the square                      % CHANGE HERE

areas = square_area_V2(alpha_P_mat(1), alpha_P_mat(2), period);   % store the areas as well
newton_steps = [0];                        % number of steps to convergence
function_eval = zeros(1,2);                       % what the function evaluates to (should be close to 0).

for i=1:1000   % gradually change the height of the rectangle
    
    heights(i+1) = 1+height_perturb*i;
    
    h = heights(i+1); % Height of parallelogram

    guesses = alpha_P_mat(i,:); % First column: init_angle and init_pos  

    for j=1:10
        guess = guesses(j,:)';  % column vector
        alpha = guess(1); P = guess(2);

        % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
        [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
        F_n = [F_alpha(period+1); F_P(period+1)];

        dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

        guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

        alpha = guess(1); P = guess(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEMP SOLUTION for when we hop out of range
        alpha = mod(alpha, pi); P = mod(P, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guesses(j+1,:) = [alpha P];

        if (norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6)  ||  (j == 10)% relative error
            
            % For function evaluation. Are we at [0,0] yet?
            [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];
            dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);
            
            if norm([F_n(2)-guesses(j+1,2) dA_dP]) < 1e-6
                areas(i+1) = parallelogram_area(h, gamma, alpha, P, period);
                newton_steps(i+1) = j;    % we converged on jth iteration
                function_eval(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
                break
            end
        end
    end
    
    alpha_P_mat(i+1,:) = guesses(end,:);  % the converged value
end

bigmat = [heights' alpha_P_mat areas' newton_steps' function_eval];
toc


%% Going BACKWARDS in height
tic

gamma = pi/2; % angle for parallelogram

period = 10; % what period are we seaching for?

eps = 1e-6;   % for finite differencing

height_perturb = 1e-3; % perturbing height of rectangle

heights_back = [1];             % The heights, filled in for square already
alpha_P_mat_back = [atan(4/1) 0.125]; % Already filled in for the square            % CHANGE HERE


areas_back = square_area_V2(alpha_P_mat_back(1), alpha_P_mat_back(2), period);   % store the areas as well
newton_steps_back = [0];                        % number of steps to convergence
function_eval_back = zeros(1,2);                       % what the function evaluates to (should be close to 0).


for i=1:800   % gradually change the height of the rectangle
    
    heights_back(i+1) = 1-height_perturb*i;
    
    h = heights_back(i+1);  % height of parallelogram

    guesses = alpha_P_mat_back(i,:); % First column: init_angle and init_pos  

    for j=1:10
        guess = guesses(j,:)';  % column vector
        alpha = guess(1); P = guess(2);

        % F^{N}(\alpha, \P) term. Where are we after a period number of bounces
        [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
        F_n = [F_alpha(period+1); F_P(period+1)];

        dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);

        guess = guess - inv(parallelogram_jacobian_area(h, gamma, alpha, P, period)) * [F_n(2)-P; dA_dP];

        alpha = guess(1); P = guess(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEMP SOLUTION for when we hop out of range
        alpha = mod(alpha, pi/2); P = mod(P, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        guesses(j+1,:) = [alpha P];

        if norm(guesses(j+1,:)- guesses(j,:))/norm(guesses(j,:)) < 1e-6  ||  (j == 10)% relative error
            
            % For function evaluation
            [F_side F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
            F_n = [F_alpha(period+1); F_P(period+1)];
            dA_dP = (parallelogram_area(h, gamma, alpha, P+eps, period)-parallelogram_area(h, gamma, alpha, P-eps, period))/(2*eps);
            
            if norm([F_n(2)-guesses(j+1,2) dA_dP]) < 1e-6
                areas_back(i+1) = parallelogram_area(h, gamma, alpha, P, period);
                newton_steps_back(i+1) = j;    % we converged on jth iteration
                function_eval_back(i+1,:) = [F_n(2)-guesses(j+1,2) dA_dP];
                break
            end
        end
    end
    
    alpha_P_mat_back(i+1,:) = guesses(end,:);  % the converged value
end

bigmat_back = [heights_back' alpha_P_mat_back areas_back' newton_steps_back' function_eval_back];
toc

%%

period10mat = [flip(bigmat_back); bigmat(2:end,:)];
%save('Data\rec_period10_alpha2_3.mat', 'period10mat');





%% This matches the result from rec_newton
plot3(period10mat(:,1), period10mat(:,2), period10mat(:,3), 'o'); hold on
set(gca,'XLim',[0 max(period10mat(:,1))], 'YLim',[0 pi/2], 'ZLim',[0 1]);
xlabel('Height of rectangle')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
yticks([0 pi/4 pi/2])
yticklabels({'0','\pi/4', '\pi/2'})
title('Continuation Period-10')

plot3(period10mat(:,1), atan(4.*period10mat(:,1)./1), period10mat(:,3), 'r-')

%% Area

plot(period10mat(:,1), period10mat(:,4), 'o')
xlabel('Height of Parallelogram')
ylabel('Area of trajectory')
title('Area vs. height of rectangle')

%% Newton steps

plot(period10mat(:,1), period10mat(:,5), 'o')
xlabel('Height of Parallelogram')
ylabel('Number of Newton steps')
title('Number of Newton steps vs. height of Parallelogram')
yticks([0 1 2 3 4 5])
yticklabels({'0','1','2','3','4','5'})
set(gca, 'YLim',[0 max(period10mat(:,5))+1]);




%% Following a period-6 orbit by keeping h constant and varying gamma.
clear all; clc
tic
period = 6*1; % what period are we seaching for?

h = 1; % keep this fixed for now

% We know for square that the angle and position that is a period-6 orbit
% and that maximises the area is angle=atan(2), pos = 0.25. This is our
% initial guess.

gamma = pi/2;             % for the square
alpha = atan(2); P = 0.25;

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
toc

% Column 1 is gamma
% Column 2 is initial angle, alpha
% Column 3 is initial position, P
% Column 4 is the area
% Column 5 is the number of Newton steps
% Column 6 & 7 are function evaluations for g^(N)(alpha,P)-P and dA/dP
% Column 8 is the function evaluation for alpha, f^(N)(alpha,P)-alpha
test_mat = test_mat(1:(end-1),:)

%%

index = size(test_mat,1);   % last row

gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

period = 36;

temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

temp_mat
        
%%
N = 50; % number of animation iterations

row = size(temp_mat,1) - 2;
h = 1; gamma = temp_mat(row, 1) ;

alpha_star = temp_mat(row, 2); P_star = temp_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%

test_mat = [test_mat; temp_mat(2:(end-1),:)];


%% The last continuation before non-convergence

N = 24; % number of animation iterations

row = size(test_mat,1) - 3;
h = 1; gamma = test_mat(row, 1) ;

alpha_star = test_mat(row, 2); P_star = test_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Animating the continuation of orbits

for j=1:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end


%% To continue the periodic orbit once it passes through the vertex. We increase the period we want to follow.
% As a period-2 orbit is also a period-6 orbit. We find that at each
% bifurcaton for this periodic orbit, the period increases by 4.

h = 1; gamma = pi/2;    % dimensions of initial parallelogram, we keep h fixed

alpha = atan(2); P = 0.25;  % the angle and position in the square for which we start following from
% (note this normally gives least period-2 orbits in the square)
    
Nmax = 36;

for i=6:6:Nmax
    
    period = i;          % following a period-i orbit
    
    if i == 6   % first iteration
        test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        test_mat = test_mat(1:(end-1),:);   % delete last row
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = test_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        test_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        
        test_mat(:,4) = test_mat(:,4)./(prod(6:4:Nmax).*h)   % normlization, to make numbers smaller
        
        
        
    else % use the previous iterations results to start from
        index = size(test_mat,1);   % last row

        gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
        alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

        temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = temp_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        temp_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        temp_mat(:,4) = temp_mat(:,4)./(prod(6:4:Nmax).*h);   % normlization, to make numbers smaller
        
        test_mat = [test_mat; temp_mat(2:(end-1),:)];  % appending the new results to the matrix

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP PLOTTING

    N = 100; % iterations for animations

    index = size(test_mat,1);
    gamma = test_mat(index, 1);
    alpha_star = test_mat(index, 2);
    P_star = test_mat(index, 3);

    [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
    P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    w = waitforbuttonpress;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


test_mat










%% 3D Plot of gamma, alpha_0, position0
plot3(test_mat(:,1), test_mat(:,2), test_mat(:,3), 'go'); hold on
set(gca,'XLim',[0 pi/2], 'YLim',[0 pi], 'ZLim',[0 2+2*h/max(test_mat(:,1))]);
xlabel('Gamma of parallelogram')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
title(sprintf('Continuation Period-%d', period))


%%

%line curve of vertice position
gammas = linspace(0.01,pi/2, 500);
plot3(gammas, 1.3.*ones(size(gammas)), 1 + 1./sin(gammas), 'k-'); hold on


%% Area

plot(test_mat(:,1), test_mat(:,4), 'o')
xlabel('gamma of parallelogram')
ylabel('Area of trajectory')
title('Area vs. gamma of parallelogram')


%% Newton steps

plot(test_mat(:,1), test_mat(:,5), 'o')
xlabel('gamma of parallelogram')
ylabel('Number of Newton steps')
title('Number of Newton steps vs. gamma of parallelogram')








%% ANOTHER PERIOD-6

clear all; clc
period = 6*4; % what period are we seaching for?

h = 1; % keep this fixed for now

gamma = pi/2;             % for the square
alpha = atan(2); P = 0.75;

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
test_mat = test_mat(1:(end-1),:)



%%

index = size(test_mat,1);   % last row

gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

period = 14;

temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

temp_mat

%test_mat = [test_mat; temp_mat(2:(end-1),:)];
        
%%
N = 50; % number of animation iterations

row = size(temp_mat,1) - 20;
h = 1; gamma = temp_mat(row, 1) ;

alpha_star = temp_mat(row, 2); P_star = temp_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%

test_mat = [test_mat; temp_mat(2:(end-1),:)];


%% The last continuation before non-convergence

N = 24; % number of animation iterations

row = size(test_mat,1) - 3;
h = 1; gamma = test_mat(row, 1) ;

alpha_star = test_mat(row, 2); P_star = test_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Animating the continuation of orbits

for j=1:20:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end











%% Period-2 orbit bifurcating into period-6 orbit.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period = 2*6*10;  % following a period-6 orbit
h = 1; gamma = pi/2;    % dimensions of parallelogram

alpha = pi/2; P = 3.5;  % the angle and position in the square for which we start following from
% (note this normally gives least period-2 orbits in the square)

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

test_mat

%% The animation for the period-2 to period-6 bifurcation

for j=5000:10:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 2,6,10,14,18,26';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end


%% To continue the periodic orbit once it passes through the vertex. We increase the period we want to follow.
% As a period-2 orbit is also a period-6 orbit. We find that at each
% bifurcaton for this periodic orbit, the period increases by 4.

h = 1; gamma = pi/2;    % dimensions of initial parallelogram, we keep h fixed

alpha = pi/2; P = 3.5;  % the angle and position in the square for which we start following from
% (note this normally gives least period-2 orbits in the square)
    
Nmax = 54;

for i=6:4:Nmax
    
    period = i;          % following a period-i orbit
    
    if i == 6   % first iteration
        test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        test_mat = test_mat(1:(end-1),:);   % delete last row
        
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = test_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        test_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        
        test_mat(:,4) = test_mat(:,4)./(prod(6:4:Nmax).*h)   % normlization, to make numbers smaller
        
        
        
    else % use the previous iterations results to start from
        index = size(test_mat,1);   % last row

        gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
        alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

        temp_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
        
        %test_mat(:,9) = test_mat(:,4);
        
        % to make all the areas comparable as we have our area func depends on period
        cum_area = temp_mat(:,4) - period .* h; % the cumulative area. Area = period*h - cum
        temp_mat(:,4) = prod(6:4:Nmax).*h - cum_area .* prod(6:4:Nmax)./period;  
        temp_mat(:,4) = temp_mat(:,4)./(prod(6:4:Nmax).*h);   % normlization, to make numbers smaller
        
        test_mat = [test_mat; temp_mat(2:(end-1),:)];  % appending the new results to the matrix

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEMP PLOTTING

    N = 100; % iterations for animations

    index = size(test_mat,1);
    gamma = test_mat(index, 1);
    alpha_star = test_mat(index, 2);
    P_star = test_mat(index, 3);

    [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
    P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off
    w = waitforbuttonpress;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


test_mat


%% Area

plot(test_mat(:,1),test_mat(:,4), 'o', 'Markersize', 3)
xlabel('\gamma'); ylabel('Area')

% Each sepearte curve is for a different period
% Eg: far most right curve is for period-2, discontinuties are for when we
% hit the vertex

%% Initial Angle

plot(test_mat(:,1),test_mat(:,2), 'o', 'Markersize', 3)
xlabel('\gamma'); ylabel('Initial Angle')

%% Initial Position

plot(test_mat(:,1),test_mat(:,3), 'o', 'Markersize', 4)
xlabel('\gamma'); ylabel('Initial Position')

%% Newton Steps

plot(test_mat(:,1),test_mat(:,5), 'o', 'Markersize', 4)
xlabel('\gamma'); ylabel('Number of Newton Steps')


%% All the different angles and positions visited 

period = Nmax;

for j = 1:(size(test_mat, 1))    % all entries are good
    gamma = test_mat(j, 1);
    alpha = test_mat(j, 2); P = test_mat(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    % adding a period column to the matrix test_mat
    next_p = find(abs(p - p(1)) < 1e-6 == 1 & abs(a - a(1)) < 1e-6 == 1, 2); % when do we begin repeating?
    test_mat(j,9) = next_p(2)-next_p(1);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end


set(gca,'XLim',[0 pi/2], 'YLim',[0 pi], 'ZLim',[0 8.5]);
xlabel('Gamma of parallelogram')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
period = prod(6:4:Nmax);
title(sprintf('Continuation Period-%d', period))

[x y] = meshgrid(0:0.01:pi/2, 0:0.01:pi); 
z1 = 1 + h./sin(x);    
z2 = 2 + h./sin(x);
z3 = 2 + 2*h./sin(x);

surf(x,y,z1,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on %Plot the surface
surf(x,y,z2,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y,z3,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, ones(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, zeros(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on

view(90, 0) % az, el


%% Plot of gamma against the period

plot(test_mat(:,1), test_mat(:,9), 'o')
xlabel('\gamma'); ylabel('Period')
set(gca,'XLim',[0 pi/2]);










%% NEW JACOBIAN TESTING. (solving the alpha-P equation (singular))
alpha_star = 1.324; P_star = 0.123; period = 4;

h = 1.236; gamma = pi/2-0.1;

parallelogram_jacobian(h, gamma, alpha_star, P_star, period)

% ONLY ONE EIGENVALUE 1 now (compared to square and rectangle)

%%
clear all; clc;

h = 1; gamma = pi/2-0.1;


init_alpha = pi/3; init_P = 0.1;
init_guess = [init_alpha; init_P]; % init_angle and init_pos

period = 4; % what period are we seaching for?

alpha = init_alpha; P = init_P;
guess=  init_guess;

for j=1:7
    disp(j)
    % F^{N}(\alpha, \P) term
   [F_alpha, F_P] = parallelogram_map(h, gamma, alpha, P, period);
    F_n = [F_alpha(period+1); F_P(period+1)];
    
    parallelogram_jacobian(h, gamma, alpha, P, period)- eye(2)
    
   guess = guess - inv(parallelogram_jacobian(h, gamma, alpha, P, period) - eye(2)) * (F_n - guess);
   
   alpha = mod(guess(1), pi); P = mod(guess(2), 1);
   
end