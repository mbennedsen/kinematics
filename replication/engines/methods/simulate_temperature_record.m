function S=simulate_temperature_record(C,d,rep)
% Generate one of the five calibrated temperature paths used in Tables S4–S6.
validateattributes(d,{'numeric'},{'scalar','integer','>=',1,'<=',numel(C.scenarios)});
assert(ismember(d,C.selected),'Design is outside the calibration.');
if ismember(d,[1 2 4])
    S=simulate_integrated_path(C,d,rep); return;
end
% Reuse exactly the observation-noise stream, with no additional state noise.
S=simulate_integrated_path(C,1,rep);
t=C.sample.t; n=numel(t); alpha=zeros(3,n);
if d==6
    P=C.parametric.quadratic; u=t-P.tRef; b=P.coeff;
    alpha(1,:)=(b(1)+b(2)*u+b(3)*u.^2)';
    alpha(2,:)=(b(2)+2*b(3)*u)';
    alpha(3,:)=2*b(3);
elseif d==7
    P=C.parametric.changepoint; u=t-P.tRef; b=P.coeff;
    alpha(1,:)=(b(1)+b(2)*u+b(3)*max(0,t-P.breakDate))';
    % Match the production estimator's left derivative at the knot.
    alpha(2,:)=(b(2)+b(3)*double(t>P.breakDate))';
    % Classical instantaneous acceleration is undefined at the velocity jump.
    alpha(3,abs(t-P.breakDate)<1e-9)=NaN;
    assert(abs(t(end)-P.breakDate)>1e-9,'Endpoint cannot coincide with the break.');
else
    error('Unknown DGP index.');
end
S.alpha=alpha;
S.y=alpha(1,:)'+P.covariateContribution+S.noise;
S.truth=[alpha(1,end),10*alpha(2,end),100*alpha(3,end), ...
    100*(alpha(2,end)-alpha(2,C.index0))/(C.endTime-C.actualT0)];
end
