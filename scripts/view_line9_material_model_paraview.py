"""Save ParaView views of coarse display grids; no solver."""
import argparse
import json
from pathlib import Path
from paraview.simple import (XMLImageDataReader, GetActiveViewOrCreate, Show, ColorBy,
    GetColorTransferFunction, GetScalarBar, Text, Sphere, Threshold, Render, SaveState, SaveScreenshot)

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--out',type=Path,required=True);p.add_argument('--dimension',choices=['2d','3d'],required=True)
p.add_argument('--revision',default='r2',choices=['r2','r3'])
args=p.parse_args();out=args.out.resolve()
state=out/f'model_{args.dimension}_{args.revision}.pvsm';screenshot=out/f'paraview_{args.dimension}_{args.revision}.png'
if state.exists() or screenshot.exists(): raise ValueError('Refuse to overwrite ParaView preview')
file='model_view_2d.vti' if args.dimension=='2d' else 'model_view_local3d.vti'
reader=XMLImageDataReader(registrationName='Line9 research material geometry',FileName=[str(out/file)])
reader.UpdatePipeline()
view=GetActiveViewOrCreate('RenderView');view.ViewSize=[1500,850]
view.Background=[.97,.97,.97];view.UseColorPaletteForBackground=0
if args.dimension=='3d':
    solid=Threshold(Input=reader);solid.Scalars=['CELLS','MaterialID'];solid.LowerThreshold=1;solid.UpperThreshold=3
    solid.UpdatePipeline();shown=solid
else: shown=reader
display=Show(shown,view);display.Representation='Surface';ColorBy(display,('CELLS','MaterialID'))
display.Ambient=1;display.Diffuse=0
lut=GetColorTransferFunction('MaterialID');lut.InterpretValuesAsCategories=1
lut.Annotations=['0','air','1','silty clay','2','mudstone','3','sandstone']
lut.IndexedColors=[217/255]*3+[244/255,177/255,131/255,191/255,144/255,0,127/255,96/255,0]
display.SetScalarBarVisibility(view,True)
bar=GetScalarBar(lut,view);bar.Title='Material';bar.ComponentTitle='';bar.TitleColor=[.1]*3;bar.LabelColor=[.1]*3
display.DataAxesGrid.GridAxesVisibility=1
display.DataAxesGrid.XTitle='Profile X (m)';display.DataAxesGrid.YTitle='Elevation (m)';display.DataAxesGrid.ZTitle='Local cross-track (m)'
for axis in ['X','Y','Z']:
    setattr(display.DataAxesGrid,axis+'TitleColor',[.15]*3)
    setattr(display.DataAxesGrid,axis+'LabelColor',[.15]*3)
view.CameraViewUp=[0,1,0];view.CameraParallelProjection=1
if args.dimension=='2d':
    view.InteractionMode='2D';view.CameraPosition=[150,442.5,500];view.CameraFocalPoint=[150,442.5,0];view.CameraParallelScale=125
    note='Complete4-material research model | X220 to25 acquisition |15m AGL midpoint\n0.1m display grid;0.025m native geometry | no FDTD execution'
else:
    view.InteractionMode='3D';view.CameraPosition=[242,481,80];view.CameraFocalPoint=[196.75,434,12];view.CameraParallelScale=40
    m=json.loads((out/'manifest.json').read_text('utf-8'));c=next(c for c in m['cases'] if c['id']=='local3d_1_cross_track')
    for key,color in [('tx_m',[.9,.15,.1]),('rx_m',[.1,.35,.85])]:
        point=c[key];s=Sphere(registrationName=key);s.Center=[point[0]+184.75,point[1]+405,point[2]];s.Radius=.35
        d=Show(s,view);d.DiffuseColor=color
    note='Local3D near ZK08: invariant extrusion 24x24m; baseline 1.3m\nDisplay grid 0.2m; native grid 0.05m; spheres mark ideal Tx/Rx\nResearch assumptions; no FDTD execution or physical acceptance'
text=Text();text.Text=note
td=Show(text,view);td.WindowLocation='Upper Left Corner';td.FontSize=11;td.Color=[.1]*3
Render(view);SaveState(str(state));SaveScreenshot(str(screenshot),view,ImageResolution=[1500,850])
print(f'ParaView {args.dimension} model state/screenshot saved',flush=True)
