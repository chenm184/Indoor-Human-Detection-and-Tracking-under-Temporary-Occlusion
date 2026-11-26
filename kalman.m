function [Xhat, P_history] = kalman(Z, detected, F, H, Q, R, x0, P0)

N = size(Z,2);
Xhat = zeros(4,N);
P_history = zeros(4,4,N);
xhat = x0;
P = P0;

for k = 1:N

    %% Prediction
    xPred = F * xhat;
    PPred = F * P * F' + Q; % estimation cova matrix

    % Update only if target is detected
    if detected(k) && ~any(isnan(Z(:,k)))

        innovation = Z(:,k) - H*xPred;
        
        K = PPred*H'*(R + H*PPred*H')^-1; % kalman gain

        xhat = xPred + K*innovation;

        P = (eye(4) - K*H)*PPred;   

    else

        % Prediction only
        xhat = xPred;
        P = PPred;

    end

    Xhat(:,k) = xhat;
    P_history(:,:,k) = P;

end

end