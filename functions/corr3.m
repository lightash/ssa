function [r12] = corr3(x1, x2)
% correlate on [ ignalsignalsigna ]
%         from [ signal---------- ]
%      through [ -----signal----- ]
%           to [ ----------signal ]

    N = length(x1);

    x1 = [x1 x1 x1];

    c = xcorr(x2, x1);

    r12 = c( N+1 : length(c) - 3*N );

end