function [Fxpts,Fypts,Txpts,Typts,varargout] = shall_not_be_denied(Fxpts,Fypts,Txpts,Typts,FprojectedBP,TprojectedBP)
%%  function [xpts,ypts] = shall_not_be_denied(xpts,ypts,projectedBP)

Fx_diff = Fxpts(1) - FprojectedBP(1);
Fy_diff = Fypts(1) - FprojectedBP(2);

Tx_diff = Txpts(1) - TprojectedBP(1);
Ty_diff = Typts(1) - TprojectedBP(2);

if (abs(Fx_diff) + abs(Fy_diff)) > (abs(Tx_diff) + abs(Ty_diff))
    Fxpts = Fxpts - Fx_diff;
    Fypts = Fypts - Fy_diff;
    varargout{1} = [Fx_diff,Fy_diff];
    varargout{2} = 'front';
else
    Txpts = Txpts - Tx_diff;
    Typts = Typts - Ty_diff;
    varargout{1} = [Tx_diff,Ty_diff];
    varargout{2} = 'top';
end

end
