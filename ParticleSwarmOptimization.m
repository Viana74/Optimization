% ############################################################################
% Subject: Paticle Swarm Optimization
% Authors: Éverton Viana, Mathias Titton
% Date: Dez./2019
% 
% Problem: Find the max value of the objective function bellow:
% f(x) = x*sin(10*pi*x)+1, where x_min <= x <= x_max com x \in R.
% ############################################################################
clc;
clear all;
format long;

np = 5;         % Number of particles belongs to swarm.
niter = 150;    % Number of iterations for each swarm.
nenx = 20;      % Number of swarms (executions).
x_min = -1;     % Minimal limit for search space.
x_max = 3;      % Maximal limit for search space.
v_min = -0.01;  % Minimal initial speed.
v_max = 0.01;   % Maximal initial speed.

% Weights w, c1, c2, r1 e r2.
k = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
X = 2*k/abs(2-phi-sqrt(phi^2-4*phi));
w = X;
c1 = X*phi1;
c2 = X*phi2;

% Animation:
Plotar = 0;         % Plot figures / Generate a GIF? 0 = No, 1 = Yes.
AnimEnx = 1;        % Define the number of the swarm (or execution) for creating the GIF.
nitergif = niter;   % Gif iterations.

% Coast function (Fitness):
fit = @(x) x.*sin(10*pi*x) + 1;

%% Swarm movement!
fprintf('Execuções: \n');
for m = 1:nenx
    % Position and initial speed for each particle of the swarm (ramdomic):
    xx(:,1,m)  = x_min + (x_max-x_min)*rand(np,1);
    vv(:,1,m)  = v_min + (v_max-v_min)*rand(np,1);

    % Best position each particle has ever passed (equal initial position):
    pp(:,1,m)  = xx(:,1,m);
    
    % Best global solution among all particles:
    [gb(:,1,m),pb] = max(fit(pp(:,1,m)));       % Best global solution.
    gg(1,1,m) = pp(pb,end,m);  

    for t = 1:niter    
        % Update the position and speed for each particle.
        for l = 1:np
            % Updating the coefficients r1 and r2.
            r1 = rand(1);
            r2 = rand(1);              
            vv(l,t+1,m) = w*vv(l,t,m)+r1*c1*(pp(l,t,m)-xx(l,t,m))+r2*c2*(gg(:,t,m)-xx(l,t,m));
            
            % Update position.
            xx(l,t+1,m) = xx(l,t,m) + vv(l,t+1,m);
        
            % Position control to avoid over limit.
            if xx(l,t+1,m) < x_min
                xx(l,t+1,m) = x_min;
            elseif xx(l,t+1,m) > x_max
                xx(l,t+1,m) = x_max;
            end
            
            % Check if the particle is better than itself previously.
            if fit(xx(l,t+1,m)) > fit(pp(l,t,m))
                pp(l,t+1,m) = xx(l,t+1,m);
            else
                pp(l,t+1,m) = pp(l,t,m);
            end
        end
        
        % Update global solution.
        [gb(:,t+1,m),pb] = max(fit(pp(:,t+1,m)));
        gg(:,t+1,m) = pp(pb,t+1,m);            
    end
    
    % Show solution.
    [best_finc(m,1),ind] = max(gb(1,:,m));
    best_xinc(m,1) = gg(1,ind,m);
    msg = [num2str(m), '  Solução: x = ',num2str(best_xinc(m,1)), '; Fitness: ',num2str(best_finc(m,1))];
    disp(msg);
end

%% Outcome:
fprintf('\n');
[~,ind] = max(best_finc);
[~,ind2] = min(best_finc);

fprintf('Best solution in the %d ª iteration.\n', ind);
fprintf('Worst solution in the %d ª iterairon.\n', ind2);
msg = ['Maximal fitness value: ',num2str(max(best_finc))];
disp(msg);
msg = ['Minimal fitness value: ',num2str(min(best_finc))];
disp(msg);  
fprintf('\n');

%% GIF generate:
if Plotar == 1
    h = figure;
    axis tight manual;   
    filename = 'animation.gif';    
    xplot = linspace(x_min,x_max,10000);
    yplot = xplot.*sin(10*pi.*xplot)+1;
    [Y,ind] = max(yplot);
    
    for j = 1:nitergif
        delete(findall(gcf,'type','annotation'));
        delete(findall(gcf,'type','text'));
        plot(xplot,yplot,'b','Linewidth',1);
        hold on;
        title('Maximization problem for f(x) = xsin(10\pix) + 1');
        grid on;
        grid minor;
        xlabel('X');
        ylabel('Fitness');
        plot(xplot(1,ind),Y(1,1),'r*');
        plot(xx(:,j,AnimEnx),fit(xx(:,j,AnimEnx)),'.k','MarkerSize',15);
        
        % Loop to avoid plot when the swarm is near of the global solution.
        limMax = Y(1,1)+0.1;
        limMin = Y(1,1)-0.1;
        for ii = 1:np
            s = [' ' num2str(ii)];
            if ~((fit(xx(ii,j,AnimEnx))<limMax)&&(fit(xx(ii,j,AnimEnx))>limMin))
                text(xx(ii,j,AnimEnx),fit(xx(ii,j,AnimEnx)),s,'Color','m','FontSize',10)
            end
        end
        legend('Function','Max Global','Swarms','Location','NorthWest');
        dim = [.4 .5 .3 .392];
        str = sprintf('Iteration: %d',j); 
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        hold off;
        
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);     
        if j == 1 
            imwrite(imind,cm,filename,'gif','DelayTime',0.75,'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','DelayTime',0.75,'WriteMode','append');
        end  
    end
end
