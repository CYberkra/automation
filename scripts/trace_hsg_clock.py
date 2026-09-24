"""Execute only official HSG choreography with symbolic clocks, no field solver."""
import ast,argparse,json,hashlib
from pathlib import Path
from types import SimpleNamespace

p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
tree=ast.parse(a.source.read_text(encoding='utf-8'))
cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CUDASubgridUpdater')
methods=[n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name in ('hsg_1','hsg_2')]
scope={};exec(compile(ast.fix_missing_locations(ast.Module(body=methods,type_ignores=[])),str(a.source),'exec'),scope)

class Trace:
    def __init__(self,r):
        self.r=r;self.iteration=0;self.n=0;self.e_level=0;self.h_level=0;self.rows=[];self.saved=[]
        self.G=SimpleNamespace();self.grid=SimpleNamespace(ratio=r,update_electric_is=self.e_is,update_magnetic_is=self.h_is,update_electric_os=lambda g:None,update_magnetic_os=lambda g:None)
        self.precursors=SimpleNamespace(update_electric=lambda:None,update_magnetic=lambda:None,
          interpolate_electric_in_time=lambda m:setattr(self,'e_level',(self.n-1)*r+m),
          interpolate_magnetic_in_time=lambda m:setattr(self,'h_level',(self.n-.5)*r+m),
          calc_exact_electric_in_time=lambda:setattr(self,'e_level',self.n*r),
          calc_exact_magnetic_in_time=lambda:setattr(self,'h_level',(self.n+.5)*r))
    def e_is(self,p):self.rows.append(dict(coarse_iteration=self.n,update='E',fine_output_index=self.iteration+1,needed_H_time=self.iteration+.5,precursor_time=self.h_level,lag=self.iteration+.5-self.h_level))
    def h_is(self,p):self.rows.append(dict(coarse_iteration=self.n,update='H',fine_E_index=self.iteration,needed_E_time=self.iteration,precursor_time=self.e_level,lag=self.iteration-self.e_level))
    def update_electric_sources(self):self.iteration+=1
    def store_outputs(self):self.saved.append(self.iteration)
    def __getattr__(self,name):
        if name.startswith('update_'):return lambda:None
        raise AttributeError(name)

# Use explicit phase-specific precursor clocks; the methods remain unmodified.
results=[]
for r in (1,3,5):
 t=Trace(r)
 for n in range(4):
  t.n=n;scope['hsg_2'](t)
  t.precursors.interpolate_magnetic_in_time=lambda m,n=n,r=r:setattr(t,'h_level',(n-.5)*r+m)
  t.precursors.calc_exact_magnetic_in_time=lambda n=n,r=r:setattr(t,'h_level',(n+.5)*r)
  t.precursors.interpolate_electric_in_time=lambda m,n=n,r=r:setattr(t,'e_level',n*r+m)
  scope['hsg_1'](t)
  # Restore hsg_2 E interpolation across previous/current main E.
  t.precursors.interpolate_electric_in_time=lambda m,r=r:setattr(t,'e_level',(t.n-1)*r+m)
  t.precursors.interpolate_magnetic_in_time=lambda m,r=r:setattr(t,'h_level',(t.n-.5)*r+m)
 assert t.saved==list(range(4*r+1)),(r,t.saved)
 assert all(row['lag']==(r-1)/2 for row in t.rows), (r,t.rows)
 results.append(dict(ratio=r,lag_in_fine_steps=(r-1)/2,lag_in_coarse_steps=(r-1)/(2*r),saved_indices=t.saved,trace=t.rows))
out=dict(solver_invoked=False,method='Unmodified official hsg_1/hsg_2 AST executed against symbolic precursor times; interpolation weights taken from official linear time blend. Initial zero-history transient is not a physical negative-time field.',source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),results=results,checks_passed=True,limitation='Clock-label consistency only, not a Maxwell solution or proof that applying a timestamp correction fixes coupling physics.')
with a.output.open('x',encoding='utf-8') as f:json.dump(out,f,indent=2)
print('Symbolic HSG clock checks passed')
