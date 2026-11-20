function bands = bootstrap_path_intervals(pointAlpha, errors, ciAlpha)
% Pointwise basic-bootstrap bands, including the joint velocity difference.
% errors is time x four estimands x replications. Never difference marginal CIs.
point = kinematic_path_values(pointAlpha);
assert(size(errors,1)==size(point,1) && size(errors,2)==4);
valid = squeeze(all(all(isfinite(errors(:,1:3,:)),1),2));
if size(errors,1)>120
    valid = valid & squeeze(all(isfinite(errors(121:end,4,:)),1));
end
valid = valid(:);
bands = struct('low',nan(size(point)), 'high',nan(size(point)), ...
    'ciLevel',1-ciAlpha, 'requestedDraws',size(errors,3), ...
    'validDraws',sum(valid), 'inferenceMode','bootstrapError');
bands.velocityErrors = reshape(errors(:,2,:),size(errors,1),size(errors,3));
bands.validPathDraws = valid;
if ~any(valid), return; end
e = errors(:,:,valid);
for k=1:4
    first = 1 + 120*(k==4);
    if first>size(point,1), continue; end
    q = prctile(e(first:end,k,:),100*[ciAlpha/2,1-ciAlpha/2],3);
    bands.low(first:end,k) = point(first:end,k)-q(:,:,2);
    bands.high(first:end,k) = point(first:end,k)-q(:,:,1);
end
end
