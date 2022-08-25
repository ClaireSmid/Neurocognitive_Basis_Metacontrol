function  LL = MB_MF_llik_P6(x,data)

% parameters
b = x(1);               % softmax inverse temperature
lr = x(2);              % learning rate
lambda = x(3);          % eligibility trace decay
w = x(4);            % mixing weight low stakes
st = x(5);              % stickiness
respst = x(6);          % stickiness

% initialization
% Qmf = ones(2,2)*0.5;
% Q2 = ones(2,1)*0.5;     % Q(s,a): state-action value function for Q-learning

% UPDATE 14/05/21: resetting these to be initialised at 4.5 rather than 0.5
Qmf = ones(2,2)*4.5;
Q2 = ones(2,1)*4.5;     % Q(s,a): state-action value function for Q-learning

Tm = cell(2,1);
Tm{1} = [0 1; 1 0];     % transition matrix
Tm{2} = [0 1; 1 0];     % transition matrix
M = [0 0; 0 0];         % last choice structure
R = [0; 0];             % last choice structure
LL = 0;

% loop through trials
for t = 1:data.N
    
    if any(data.timeout(t,:))
        continue
    end
    
    if (data.stimuli(t,1) == 2) || (data.stimuli(t,1) == 4)
        R = flipud(R);
    end
    
    s1 = data.s(t,1);
    s2 = data.s(t,2);
    a = data.choice(t);
    action = a;
    a = a - (s1 == 2)*(2);
    
    Qmb = Tm{s1}'*Q2;                                           % compute model-based value function
    
%     if data.stake(t) == 1
%         w = w_lo;
%     else
%         w = w_hi;
%     end
    
    Q = w*Qmb + (1-w)*Qmf(s1,:)' + st.*M(s1,:)' + respst.*R;    % mix TD and model value
    
    LL = LL + b*Q(a) - logsumexp(b*Q);
    
    M = zeros(2,2);
    M(s1,a) = 1;                                                % make the last choice sticky
    
    R = zeros(2,1);
    if action == data.stimuli(t,1)
        R(1) = 1;                                               % make the last choice sticky
    else
        R(2) = 1;
    end
    
    dtQ(1) = Q2(s2) - Qmf(s1,a);                                % backup with actual choice (i.e., sarsa)
    Qmf(s1,a) = Qmf(s1,a) + lr*dtQ(1);                          % update TD value function
    
    dtQ(2) = data.points(t) - Q2(s2);                           % prediction error (2nd choice)
    
    Q2(s2) = Q2(s2) + lr*dtQ(2);                                % update TD value function
    Qmf(s1,a) = Qmf(s1,a) + lambda*lr*dtQ(2);                   % eligibility trace
    
end