
function [detected,Y]= np_detector(Ak,sigman)
N = length(Ak);
wk = sigman * randn(N,1);
Y = Ak + wk;
pfa = 0.05;
snrThdB = npwgnthresh(pfa,1,'real');
snrThNP = sqrt(db2pow(snrThdB));
gamma = snrThNP*sigman; % times normolization
detected = Y > gamma; % Y is the test statistics T(Y) > gamma ?

end
