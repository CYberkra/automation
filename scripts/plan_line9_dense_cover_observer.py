"""Geometry-only plan: interpolate observer locations, never geology or solver fields."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from line9_cover_angular_transport import cover_permittivity,C0


def main(a):
    m=json.loads((a.package/'manifest.json').read_text('utf-8'));assert m['receiver_count']==1231
    with h5py.File(a.package/m['groups'][0]['geometry']) as h:geo=h['data'][:,:,0]
    existing={tuple(p['anchors'][0]['coord']) for p in m['probes']};new=[];rx=1232
    for y in m['cover_plane_y_m']:
        j=round(y/.025)
        for x in np.arange(154,174.001,.125):
            i=round(x/.025)
            if (i,j,0) in existing:continue
            assert geo[i-1:i+2,j-1:j+2].shape==(3,3) and np.all(geo[i-1:i+2,j-1:j+2]==1)
            anchors=[]
            for k,(coord,outputs) in enumerate([([i,j,0],['Ez','Hx','Hy']),([i,j-1,0],['Hx']),([i-1,j,0],['Hy'])]):
                anchors.append({'receiver_index':rx,'coord':coord,'position_m':[coord[0]*.025,coord[1]*.025,0.],
                                'name':f'p_dense_cover_j{j}_i{i}_a{k}','outputs':outputs});rx+=1
            new.append({'id':f'dense_cover_j{j}_i{i}','x_m':float(x),'y_m':j*.025,'anchors':anchors})
    assert len(new)==360 and rx==2312
    f=np.linspace(20e6,170e6,501);k=2*np.pi*f/C0*np.sqrt(cover_permittivity(f))
    proposal={'status':'DENSE_COVER_GEOMETRY_VERIFIED_PROPOSAL_NOT_PREPARED_FROZEN_OR_RUN',
              'parent_manifest_sha256':sha(a.package/'manifest.json'),'generator_sha256':sha(__file__),
              'parent_native_sha256':'a601791ce1cd85585d9f5cad24fd026ccdd03baa3e55cf1b0ff7f281cec3ac26',
              'new_probes':new,'preserve_logical_points':410,'new_logical_points':360,'total_logical_points':770,
              'preserve_raw_receivers':1231,'new_raw_receivers':1080,'total_raw_receivers':2311,
              'cover_plane_y_m':m['cover_plane_y_m'],'x_range_m':[154,174],'target_observer_spacing_m':.125,
              'solver_spacing_m_unchanged':.025,'source_geometry_material_previous_channels_invariant_required':True,
              'device_history_bytes':2311*20352*6*8,'saved_history_bytes':20352*8*(1+5*770),
              'conservative_device_estimate_bytes':8400*1700*420+2*2**30+2311*20352*6*8,
              'cover_max_propagating_real_k_per_m':float(k.real.max()),'old_observer_Nyquist_per_m':float(np.pi/.5),
              'new_observer_Nyquist_per_m':float(np.pi/.125),'all_new_collocation_support_uniform_cover':True,
              'scientific_question':'Separate denser sampled full-cover propagation from restricted air-sector transport; compare retained .5m positions and explicit decimation; inspect high-q/alias versus finite aperture and heterogeneous exterior.',
              'analysis_requirements':['Full complex Debye loss counted once','All previous2051 raw receiver channels and source bitwise invariant',
                  'Full501-tone Hann and Blackman mainRx invariant','Uniform-cover sector incl. grazing safeguards and evanescent limits',
                  'Spatial aperture sensitivity and explicit .5m decimation','Independent sum checks and analytic controls','Chinese raw and split plots with common scales'],
              'limits':'Local plane support verified only; evanescent/finite aperture/exterior effects remain. No automatic solve or whole-line restart.'}
    assert not a.out.exists();a.out.write_text(json.dumps(proposal,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:proposal[k] for k in ['status','new_logical_points','total_raw_receivers','cover_max_propagating_real_k_per_m','old_observer_Nyquist_per_m','new_observer_Nyquist_per_m','conservative_device_estimate_bytes']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
