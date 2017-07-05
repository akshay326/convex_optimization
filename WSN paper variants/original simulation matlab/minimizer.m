delta = 0.8;
N = 2; % No of players
l = 2; % No of actions
W = 30; % No of states

% For now, I'll start with single player
cvx_begin
    variables sig(W) v(W);
    PI = zeros(W,W);
    sig_joint = sig(w)*(1-sig(w));
    
    for i = 1:W
        for j = 1:W
            switch i-j
                case -1
                    PI(i,j) = j*sig_joint/W;
                case 0
                    PI(i,j) = j/W + (1-2*(j/W))*sig_joint;
                case 1
                    PI(i,j) = (1-(j/W))*(1-sig_joint);
                otherwise
                    PI(i,j) = 0;
            end
        end
    end
    
    for w = 1:W
        
        minimize(norm(v*PI,Inf)) %xi1(w,sig_joint)
        subject to
           0 <= sig <= 1;
    end  
cvx_end