function [compatible,differences]=compare_fit_settings(requested,saved)
% Only compare completed whole-fit caches. Never use this to merge partial runs.
% Moving output or requesting other workers does not change a stored result.
% Keep the original key and fit.options unchanged as provenance of its draws.
requested=scientific_key(requested);saved=scientific_key(saved);
differences=compare_fields(requested,saved,'key');compatible=isempty(differences);
end

function key=scientific_key(key)
if ~isstruct(key)||~isscalar(key)||~isfield(key,'options'),return;end
if isfield(key.options,'outputDir'),key.options=rmfield(key.options,'outputDir');end
if isfield(key.options,'bootstrap')
 for f={'useParallel','nWorkers'}
  if isfield(key.options.bootstrap,f{1}),key.options.bootstrap=rmfield(key.options.bootstrap,f{1});end
 end
end
end

function out=compare_fields(a,b,path)
out={};
if isstruct(a)&&isscalar(a)&&isstruct(b)&&isscalar(b)
 fields=union(fieldnames(a),fieldnames(b));
 for i=1:numel(fields)
  f=fields{i};p=[path '.' f];
  if ~isfield(a,f),out{end+1}=sprintf('%s: present only in saved key',p); %#ok<AGROW>
  elseif ~isfield(b,f),out{end+1}=sprintf('%s: present only in requested key',p); %#ok<AGROW>
  else,out=[out,compare_fields(a.(f),b.(f),p)]; %#ok<AGROW>
  end
 end
elseif ~isequaln(a,b)
 out={sprintf('%s: saved=%s; requested=%s',path,brief(b),brief(a))};
end
end

function s=brief(value)
try,s=jsonencode(value);catch,s=['<' class(value) '>'];end
if numel(s)>180,s=[s(1:177) '...'];end
end
