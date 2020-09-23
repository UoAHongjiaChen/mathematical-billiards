h = 1;
gamma = pi/2;
period = 6;

init_alpha_1 = atan(2);

bigmat_0p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 0.25, period);
bigmat_0p75 = parallelogram_newton_solver(h, gamma, init_alpha_1, 0.75, period);
bigmat_2p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 2.25, period);
bigmat_2p75 = parallelogram_newton_solver(h, gamma, init_alpha_1, 2.75, period);

% bigmat_1p25 = parallelogram_newton_solver(h, gamma, init_alpha_1, 1.25, period);

init_alpha_2 = atan(1/2);

bigmat_1p5 = parallelogram_newton_solver(h, gamma, init_alpha_2, 1.5, period);
bigmat_3p5 = parallelogram_newton_solver(h, gamma, init_alpha_2, 3.5, period);


%%

for j = 1:(size(bigmat_0p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p25(j, 1);
    alpha = bigmat_0p25(j, 2); P = bigmat_0p25(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, p, 'bo', 'MarkerSize', 1); hold on
    %plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'bo', 'MarkerSize', 1); hold on
end

for j = 1:(size(bigmat_0p75, 1)-1)    % last entry is the one we terminate at, so is scuffed
    gamma = bigmat_0p75(j, 1);
    alpha = bigmat_0p75(j, 2); P = bigmat_0p75(j, 3);
    
    [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
    
    plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'go', 'MarkerSize', 2); hold on
end

% for j = 1:(size(bigmat_2p25, 1)-1)    % last entry is the one we terminate at, so is scuffed
%     gamma = bigmat_2p25(j, 1);
%     alpha = bigmat_2p25(j, 2); P = bigmat_2p25(j, 3);
%     
%     [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
%     
%     plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'c*', 'MarkerSize', 1); hold on
% end

% for j = 1:(size(bigmat_2p75, 1)-1)    % last entry is the one we terminate at, so is scuffed
%     gamma = bigmat_2p75(j, 1);
%     alpha = bigmat_2p75(j, 2); P = bigmat_2p75(j, 3);
%     
%     [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
%     
%     plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'm+', 'MarkerSize', 1); hold on
% end

% for j = 1:(size(bigmat_3p5, 1)-1)    % last entry is the one we terminate at, so is scuffed
%     gamma = bigmat_3p5(j, 1);
%     alpha = bigmat_3p5(j, 2); P = bigmat_3p5(j, 3);
%     
%     [s, a, p] = parallelogram_map(h, gamma, alpha, P, period);
%     
%     plot3(gamma * ones(size(a)), a, mod(p,2+2*h/sin(gamma)), 'ro', 'MarkerSize', 1); hold on
% end




%%

set(gca,'XLim',[1.35 pi/2], 'YLim',[0 pi], 'ZLim',[0 4]);
xlabel('Gamma of parallelogram')
ylabel('\alpha_0 - angle')
zlabel('P_0 - position')
title(sprintf('Continuation Period-%d', period))

%%

%gammas = linspace(0.01,pi/2, 500);
%plot3(gammas, 1.3.*ones(size(gammas)), 1 + 1./sin(gammas), 'k-'); hold on


[x y] = meshgrid(0:0.01:pi/2, 0:0.01:pi); 
z1 = 1 + h./sin(x);    
z2 = 2 + h./sin(x);
z3 = 2 + 2*h./sin(x);

surf(x,y,z1,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on %Plot the surface
surf(x,y,z2,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y,z3,'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, ones(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on
surf(x,y, zeros(size(x)),'FaceAlpha',0.5, 'Linestyle', 'None', 'FaceColor', 'r'); hold on





view(0, 0) % az, el


%%

[s, a, p] = parallelogram_map(1, pi/2, atan(2), 0.25, period);
plot3(pi/2 * ones(size(a)), a, p, 'ro', 'MarkerSize', 5); hold on





%% TESTING DELETE











%% TRYING TO FOLLOW SOME MORE

j = (size(bigmat_0p25, 1)-1);
gamma = bigmat_0p25(j, 1);
alpha = bigmat_0p25(j, 2); P = bigmat_0p25(j, 3);
[s, a, p] = parallelogram_map(h, gamma, alpha, P, period);


index = 5;

test_mat = parallelogram_newton_solver(h, gamma, a(index)+0.01, p(index)+0.01, period);

test_mat




%% FOLLOWING REPEATED Period-6
row = 850;
bigmat(row,:)

%% Testing, feel free to change this
N = 8;

h = 1; gamma = bigmat(row, 1) ; 
alpha = bigmat(row, 2); P = bigmat(row, 3);  

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

test_mat



%%

index = 200; % size(test_mat,1)-10;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);















%% Following randoms

period = 96;
h = 1; gamma = pi/2; %bigmat(row, 1) ; 
alpha = pi/2; P = 3.5;  

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);

test_mat

%%
N = 100; % iterations for animations

index = size(test_mat,1)-1;
gamma = test_mat(index, 1);
alpha_star = test_mat(index, 2);
P_star = test_mat(index, 3);

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star-0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);

%%
index = size(test_mat,1)-1

alpha_star = test_mat(index, 2); P_star = test_mat(index, 3)-0.4;

[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);




%%
























%% TESTING OTHER PERIOD-2





%% ANOTHER PERIOD-6

clear all; clc
period = 2*12; % what period are we seaching for?

h = 1; % keep this fixed for now

gamma = pi/2;             % for the square
alpha = pi/2; P = 1.5;

test_mat = parallelogram_newton_solver(h, gamma, alpha, P, period);
test_mat = test_mat(1:(end-1),:)



%%

index = size(test_mat,1);   % last row

gamma = test_mat(index,1);    % the last gamma that we converged, this is where we start from now
alpha = test_mat(index,2); P = test_mat(index,3);  % the angle and position in for which we last converged

period = 18;

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

row = size(test_mat,1) - 10;
h = 1; gamma = test_mat(row, 1) ;

alpha_star = test_mat(row, 2); P_star = test_mat(row, 3);   
[side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star+0.00, N);
P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position);



%% Animating the continuation of orbits

for j=600:10:size(test_mat,1)  % iterating over each row in the matrix
    if ~isnan(test_mat(j, 2))
        N = 50; % iterations for animations

        h = 1; gamma = test_mat(j, 1);

        alpha_star = test_mat(j, 2); P_star = test_mat(j, 3);   

        [side alpha position] = parallelogram_map(h, gamma, alpha_star, P_star, N);
        
        % adding a period column to the matrix test_mat
        next_p = find(abs(position - position(1)) < 1e-6 == 1 & abs(alpha - alpha(1)) < 1e-6 == 1, 2); % when do we begin repeating?
        test_mat(j,9) = next_p(2)-next_p(1);
    
        P_anim = parallelogram_billiards_draw(h, gamma, side, alpha, position); hold off

        Image = getframe(gcf);
%         folder = 'C:/Users/Work/Desktop/UoA/Year 4/Honours Project/Mathematical Billiards/MATLAB Code/Images/Parallelogram/Animation/Period 6';
%         baseFileName = sprintf('test%d.jpg', j);
%         fullFileName = fullfile(folder, baseFileName);
%         imwrite(Image.cdata, fullFileName);
    end
end




