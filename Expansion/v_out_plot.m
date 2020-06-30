K = 20;
v_outs = zeros(K,1);


for d = 1:K
    floor((d/K)*2^6)
    v_out = beam_ball_MOL_func(floor((d/K)*2^6));
    v_outs(d) = v_out;
end

scatter((1:K)/K,v_outs/10,'LineWidth',3)
ylim([0,1])
xlabel('Impact Point on Racket')
ylabel('ACOR')
title('ACOR no-hysteresis, one clamp right free')