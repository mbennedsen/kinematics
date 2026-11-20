function check_smooth_intervals(R,S)
for m=1:4
 if ~R.success(m),continue;end
 keep=1:4;if m==4,keep=[1 2 4];assert(R.point(m,3)==0 && isnan(R.low(m,3)) && isnan(R.high(m,3)));end
 assert(all(isfinite([R.point(m,:),R.low(m,keep),R.high(m,keep)])));
 assert(all(R.low(m,keep)<=R.high(m,keep)));
 if m==3,continue;end
 d=R.draws{m};assert(isequal(size(d),[S.B 4]));
 for k=keep
  good=d(isfinite(d(:,k)),k);assert(numel(good)==R.validDraws(m,k) && numel(good)>=ceil(.95*S.B));
  if m==2
   values=sort(good(:));pos=1+(numel(values)-1)*[.05 .95];lo=floor(pos);hi=ceil(pos);
   q=values(lo)'+(pos-lo).*(values(hi)'-values(lo)');
  else
   q=prctile(good,[5 95]);
  end
  if m==1,expected=R.point(m,k)-q([2 1]);else,expected=q;end
  assert(max(abs(expected-[R.low(m,k) R.high(m,k)]))<1e-9,'Interval reconstruction failed.');
 end
end
end
