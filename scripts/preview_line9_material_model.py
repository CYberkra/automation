"""Make physical-coordinate previews of a prepared private model; no solver."""
import argparse
import json
from pathlib import Path

import numpy as np

from build_pdf_profile_geometry import digest, save_json
from prepare_line9_material_model import classify


def main(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
    targets=[out/'model_acquisition.png',out/'model_view_2d.vti',out/'model_view_local3d.vti',out/'preview_validation.json']
    if any(p.exists() for p in targets):
        raise ValueError('Refuse to overwrite model previews')
    m=json.loads((out/'manifest.json').read_text('utf-8'))
    src=json.loads((out/'contacts_m.json').read_text('utf-8'))
    contacts=src['points_x_elevation_m']
    records=[]
    for name,stride,target in [('full2d',4,targets[1]),('local3d_1',4,targets[2])]:
        g=m['geometries'][name]
        native_dl=g['grid_xyz_m'][0]
        dl=native_dl*stride
        native_shape=np.array(g['shape_nxyz'])
        shape=native_shape//stride
        shape[2]=max(1,shape[2])
        spacing=[dl,dl,native_dl if native_shape[2]==1 else dl]
        xmin=g['profile_x_range_m'][0]
        ymin=g['elevation_range_m'][0]
        zmin=g['cross_track_range_m'][0]
        x=xmin+(np.arange(shape[0])+.5)*dl
        y=ymin+(np.arange(shape[1])+.5)*dl
        plane=classify(x,y,contacts)
        volume=np.repeat(plane[:,:,None],int(shape[2]),axis=2)
        grid=vtk.vtkImageData()
        grid.SetDimensions(*(shape+1))
        grid.SetSpacing(*spacing)
        grid.SetOrigin(xmin,ymin,zmin)
        array=numpy_to_vtk(volume.ravel(order='F'),deep=True,array_type=vtk.VTK_SHORT)
        array.SetName('MaterialID')
        grid.GetCellData().SetScalars(array)
        w=vtk.vtkXMLImageDataWriter(); w.SetFileName(str(target)); w.SetInputData(grid)
        if w.Write()!=1: raise ValueError('Preview VTK write failed')
        reader=vtk.vtkXMLImageDataReader(); reader.SetFileName(str(target)); reader.Update()
        if not np.array_equal(vtk_to_numpy(reader.GetOutput().GetCellData().GetScalars()),volume.ravel(order='F')):
            raise ValueError('Preview VTK order round-trip failed')
        records.append(dict(file=target.name,sha256=digest(target),display_spacing_xyz_m=spacing,
                            native_geometry_sha256=g['sha256'],bounds=reader.GetOutput().GetBounds(),
                            scope='Coarser display grid evaluated from the same reviewed contacts; not solver voxel grid'))
        if name=='full2d': fullplane=plane
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei']
    plt.rcParams['axes.unicode_minus']=False
    colors=['#d9d9d9','#f4b183','#bf9000','#7f6000']
    fig,axes=plt.subplots(2,1,figsize=(18,9),layout='constrained')
    full=m['geometries']['full2d']
    extent=[*full['profile_x_range_m'],*full['elevation_range_m']]
    native_cases=[c for c in m['cases'] if c['id'].startswith('full2d_s')]
    flight_x=np.array([c['profile_x_m'] for c in native_cases])
    flight_y=np.array([c['tx_m'][1]+full['elevation_range_m'][0] for c in native_cases])
    for ax in axes:
        ax.imshow(fullplane.T,origin='lower',extent=extent,aspect='equal',interpolation='nearest',
                  cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(-.5,4.5),4))
        for name,points in contacts.items():
            points=np.asarray(points); ax.plot(points[:,0],points[:,1],color='#504b43',lw=.7)
        ax.plot(flight_x,flight_y,color='#1769aa',lw=1.8,label='收发中点轨迹：离地约15m')
        for px,label in [(220,'起点 X220'),(25,'终点 X25')]:
            iy=int(np.argmin(abs(flight_x-px)))
            ax.scatter([px],[flight_y[iy]],c=['#1769aa'],s=30)
            ax.annotate(label,(px,flight_y[iy]),xytext=(0,12),textcoords='offset points',ha='center',fontsize=10)
        ax.set(ylabel='高程 (m)',ylim=full['elevation_range_m'])
        ax.grid(alpha=.15)
    axes[0].set(xlim=full['profile_x_range_m'],xlabel='原剖面横坐标 X (m)',
                title='完整研究模型：PDF岩性接触线 + 四类电性；横纵等比例')
    axes[0].axvspan(-50,0,color='#dddddd',alpha=.25)
    axes[0].text(-25,414,'端点外延',ha='center',fontsize=10)
    axes[0].legend(loc='upper right')
    axes[1].set(xlim=(220,25),xlabel='原剖面横坐标 X (m)，向右为实际采集顺序',
                title='B-scan测线区域：X220 → X25，s=220−X；0.5m站距、391站')
    fig.legend(handles=[Patch(facecolor=c,label=n) for c,n in zip(colors,['空气','粉质黏土','泥岩','砂岩'])],
               loc='outside lower center',ncol=4)
    fig.savefig(targets[0],dpi=140); plt.close(fig)
    save_json(targets[3],dict(status='PASS_PHYSICAL_COORDINATE_PREVIEW',calls_solver=False,
        manifest_sha256=digest(out/'manifest.json'),script_sha256=digest(Path(__file__)),
        png_sha256=digest(targets[0]),vtk=records,vertical_exaggeration=1,
        source_midpoint_trajectory=True,antenna_model=False,
        note='Upper view X increases; lower view follows acquisition X decreases.3D preview is local invariant extrusion near ZK08. Display meshes are coarser than native solver geometries.'))
    print('Equal-scale model/acquisition PNG and two VTK previews PASS')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args().out)
