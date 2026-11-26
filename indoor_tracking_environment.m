clear;
clc;
close all;

%% Define Variables

roomWidth  = 6;   % meters
roomHeight = 6;   % meters

% Fixed sensor position
sensorPos = [0.0 , 3.0];  
% random obstacle
obsWidth = 0.6;
obsHeight = 1.2;
dt = 0.1;              % sampling time
T  = 15;               % total simulation time
t  = 0:dt:T;
F = [1 0 dt 0;         %state model
     0 1 0 dt;
     0 0 1  0;
     0 0 0  1];
N = length(t);
sigma_a = 0.15;        %variace of motion moise
G = [0.5*dt^2  0;
     0           0.5*dt^2;
     dt          0;
     0           dt];

H = [1 0 0 0;
     0 1 0 0];        

R = diag([0.04 0.04]);   %  sensor variance = 0.04 -> std = 0.2 m

Z = nan(2,N);

%% Define Human Trajectory

% Start and end positions
startPos = [1.0, 1.0];
%endPos   = [5.0, 5.0];
humanRadius = 0.25;

% Initialize State Space
StateTrue = zeros(4,N);
StateTrue(:,1) = [startPos(1),startPos(2),0.33,0.33];


for k = 1:N-1
    
    ak = sigma_a * randn(2,1);
    StateTrue(:,k+1) = F*StateTrue(:,k) + G*ak;
end
%% Random Obstacle

obstacleCreated = false;
for k = 1:N

    if t(k) == 5 && ~obstacleCreated

        targetPos = StateTrue(1:2,k)';

        midPoint = (sensorPos + targetPos)/2;
        obstacle = [midPoint(1)-obsWidth/2, ...
                    midPoint(1)+obsWidth/2, ...
                    midPoint(2)-obsHeight/2, ...
                    midPoint(2)+obsHeight/2];

        obstacleCreated = true;
    end

end


% check whether the target is visible
A0 = 2; % attenuation factor
occlusionRatio = zeros(N,1);
isVisible = true(N,1);
numRays = 20;
for k = 1:N
    targetCenter = StateTrue(1:2,k)';
    blockedCount = 0;
    angles = linspace(0,2*pi,numRays);
    for j = 1:numRays
        samplePoint = targetCenter + humanRadius*[cos(angles(j)), sin(angles(j))];

        if lineIntersectsRectangle(sensorPos, samplePoint, obstacle)

            blockedCount = blockedCount + 1;
        end
    end
    occlusionRatio(k) = blockedCount / numRays;
end

%% Appling NP detector
attenuation = 1 - occlusionRatio;
Ak = A0 .* attenuation;
[detected,Y] = np_detector(Ak,R(1));

targetPresent = Ak > 0;
trueDetection = detected & targetPresent;
missedDetection = ~detected & targetPresent;
falseAlarm = detected & ~targetPresent;
correctRejection = ~detected & ~targetPresent;
faIndex = find(falseAlarm);
faTime = t(faIndex);

for k = 1:N
    
    if trueDetection(k)

        vk = chol(R,'lower') * randn(2,1);
        Z(:,k) = H*StateTrue(:,k) + vk;
    elseif falseAlarm(k)
        cx = roomWidth*rand();
        cy = roomHeight*rand();
        Z(:,k) = [cx;cy];
       
    end

end


%% Appling KF


Q = sigma_a^2 * (G*G'); % covariance of dynamic model 4x4
x0 = [startPos(1),startPos(2),0.0,0.0]';
P0 = diag([0.05 0.05 0.5 0.5]);
[Xhat, P_history] = kalman(Z, detected, F, H, Q, R, x0, P0);


%% Plot Static Environment

%Update xTrue, yTrue for plot
xTrue = StateTrue(1,:);
yTrue = StateTrue(2,:);

figure;
hold on;
grid on;
axis equal;

xlim([0 roomWidth]);
ylim([0 roomHeight]);

xlabel('X Position [m]');
ylabel('Y Position [m]');
title('2D Indoor Human Tracking Environment');

% Room boundary
rectangle('Position',[0 0 roomWidth roomHeight], ...
          'EdgeColor','k', ...
          'LineWidth',2);

% Obstacle
rectangle('Position',[obstacle(1), ...
                      obstacle(3), ...
                      obstacle(2)-obstacle(1), ...
                      obstacle(4)-obstacle(3)], ...
          'FaceColor',[0.5 0.5 0.5], ...
          'EdgeColor','k');

% Sensor
plot(sensorPos(1),sensorPos(2), ...
     'bs', ...
     'MarkerSize',10, ...
     'MarkerFaceColor','b');

text(sensorPos(1)+0.1,sensorPos(2)+0.15,'Sensor');

% True trajectory
plot(xTrue,yTrue,'k--','LineWidth',1.2);

legend('Room','Obstacle','Sensor','True Trajectory');

%% 6. Animate Human Motion

humanMarker = plot(xTrue(1),yTrue(1), ...
                   'ro', ...
                   'MarkerSize',10, ...
                   'MarkerFaceColor','r');
measurementMarker = plot(nan,nan, ...
                         'bo', ...
                         'MarkerSize',8, ...
                         'MarkerFaceColor','b');



losLine = plot([sensorPos(1) xTrue(1)], ...
               [sensorPos(2) yTrue(1)], ...
               'g-', ...
               'LineWidth',1.5);

statusText = text(3.8,5.5,'','FontSize',12);

for k = 1:N

    % Update human position
    humanMarker.XData = xTrue(k);
    humanMarker.YData = yTrue(k);

    %Update sensor measurement position
     measurementMarker.XData = Z(1,k);
     measurementMarker.YData = Z(2,k);

    % Update sensor-target LOS line
    losLine.XData = [sensorPos(1), Z(1,k)];
    losLine.YData = [sensorPos(2), Z(2,k)];

    % Change display depending on visibility
    if isVisible(k)

        humanMarker.MarkerFaceColor = 'r';
        losLine.Color = 'g';

        statusText.String = 'Target Visible';

    else

        humanMarker.MarkerFaceColor = 'y';
        losLine.Color = 'r';

        statusText.String = 'Target Occluded';

    end

    drawnow;

    pause(0.03);

end




%% Plotting the results
% p1:2D trajectory comparison
figure;
hold on;
grid on;
axis equal;

plot(StateTrue(1,:), StateTrue(2,:),'k--','LineWidth',1.5);

scatter(Z(1,:), Z(2,:),'filled');

plot(Xhat(1,:), Xhat(2,:),'LineWidth',1.5);

xlabel('X Position [m]');
ylabel('Y Position [m]');
title('True Trajectory, Measurements, and Kalman Estimate');

legend('True Trajectory', 'Sensor Measurements', 'Kalman Estimate');


% p2: position error
position_error = sqrt((Xhat(1,:) - StateTrue(1,:)).^2 + (Xhat(2,:) - StateTrue(2,:)).^2);
figure;
hold on;
grid on;

plot(t,position_error,'LineWidth',1.5);
xlabel('Time [s]');
ylabel('Position Error [m]');
ylim([0,1.0]);

title('Kalman Tracking Position Error');



% p3: Kalman Position Uncertainty
%{
positionUncertainty = zeros(1,N);

for k = 1:N
    positionUncertainty(k) = ...
        trace(P_history(1:2,1:2,k));
end

figure;

plot(t, positionUncertainty, 'LineWidth', 1.5);

grid on;
xlabel('Time [s]');
ylabel('tr(P_{position})');
title('Kalman Position Uncertainty');
%}
%% -------------------------------------------------------
% Local Function
% Check whether a line segment intersects a rectangle
%% -------------------------------------------------------

function intersects = lineIntersectsRectangle(p1,p2,rect)

xmin = rect(1);
xmax = rect(2);
ymin = rect(3);
ymax = rect(4);

% Rectangle edges
edges = [
    xmin ymin xmax ymin;
    xmax ymin xmax ymax;
    xmax ymax xmin ymax;
    xmin ymax xmin ymin
];

intersects = false;

for i = 1:4

    q1 = edges(i,1:2);
    q2 = edges(i,3:4);

    if segmentsIntersect(p1,p2,q1,q2)

        intersects = true;
        return;

    end

end

end


function result = segmentsIntersect(p1,p2,q1,q2)

d1 = direction(q1,q2,p1);
d2 = direction(q1,q2,p2);
d3 = direction(p1,p2,q1);
d4 = direction(p1,p2,q2);

result = (d1*d2 < 0) && (d3*d4 < 0);
end


function d = direction(a,b,c)

d = (c(1)-a(1))*(b(2)-a(2)) - (c(2)-a(2))*(b(1)-a(1));

end