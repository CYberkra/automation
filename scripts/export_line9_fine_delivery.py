"""Export approved compact delivery; leave the 90 GiB native fields in place."""
import argparse
import json
from pathlib import Path
import zipfile
from hs_capsule_identity import sha256 as sha


def read(p):return json.loads(p.read_text('utf-8'))


def export(root):
    c=read(root/'execution/execution_contract.json');v=read(root/'execution/completed_verification.json')
    assert v['completed'] and len(v['groups'])==6
    sources={}
    for name in ['execution_contract.json','completed_verification.json','execution.jsonl','preflight_verification.json','input_verification.json']:
        sources[name]=root/'execution'/name
    sources['preparation_manifest.json']=root/'package/manifest.json'
    for g in c['groups']:
        folder=Path(g['input']).parent;id=g['id'];ov=read(folder/'observer_verification.json')
        for name in ['profile.h5','observer_verification.json','snapshot_storage.json','stdout.log','stderr.log','profile.in']:
            sources[f'native/{id}/{name}']=folder/name
        targets=[150,200,240,280,320,360,400,440] if g['kind']=='local' else [100,200,400,600,1200,1800,2300]
        for t in targets:
            row=min(ov['snapshots'],key=lambda r:abs(r['iteration']*g['dt_s']*1e9-t))
            sources[f"keys/{id}/{row['file']}"]=folder/'profile_snaps'/row['file']
        sources[f'analysis_arrays/{id}_arrays.h5']=root/'analysis_arrays_r1'/(id+'_arrays.h5')
    private=root/'compact_native_retrieval_r1.zip';public=root/'compact_public_delivery_r1.zip'
    assert not private.exists() and not public.exists()
    rows=[]
    with zipfile.ZipFile(private,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for rel,p in sources.items():
            assert p.is_file();rows.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p)));z.write(p,rel)
        z.writestr('retrieval_manifest.json',json.dumps(dict(files=rows,native_snapshots_retained_on_ROG=True),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(archive=private.name,bytes=private.stat().st_size,sha256=sha(private),files=len(rows))),flush=True)
    files=[];excluded=[];reduced=[]
    compact=root/'compact_previews_r1';assert not compact.exists();compact.mkdir()
    for p in sorted((root/'delivery_r1').iterdir()):
        if p.suffix=='.gif' and p.stem.endswith('_preview') and p.stat().st_size>25*2**20:
            from PIL import Image,ImageDraw,ImageFont,ImageSequence
            images=[]
            with Image.open(p) as source:
                duration=source.info['duration']
                for im in ImageSequence.Iterator(source):
                    rgb=im.convert('RGB');crop=rgb.crop((0,35,rgb.width,rgb.height))
                    small=crop.resize((crop.width*2//3,crop.height*2//3),Image.Resampling.BOX)
                    view=Image.new('RGB',(small.width,small.height+32),'white');view.paste(small,(0,32))
                    ImageDraw.Draw(view).text((5,4),'紧凑预览：原生图像约缩小3倍；沿用预览帧，完整原生动画在ROG',fill='black',font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',14))
                    images.append(view);rgb.close();crop.close();small.close()
            dest=compact/p.name
            images[0].save(dest,save_all=True,append_images=images[1:],duration=duration,loop=0,optimize=False)
            for im in images:im.close()
            assert dest.stat().st_size<=40*2**20
            reduced.append(dict(file=dest.name,source_sha256=sha(p),source_remote_path=str(p),sha256=sha(dest),
                bytes=dest.stat().st_size,frames=len(images),duration_ms=duration,total_spatial_display_decimation_approx=3,
                temporal_sequence_unchanged_from_preview=True,display_only=True))
            files.append(dest)
        elif p.suffix=='.gif' and p.stat().st_size>40*2**20:
            excluded.append(dict(file=p.name,bytes=p.stat().st_size,sha256=sha(p),remote_path=str(p)))
        else:files.append(p)
    rows=[dict(path=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in files]
    with zipfile.ZipFile(public,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for p in files:z.write(p,p.name)
        z.writestr('preview_reduction_receipts.json',json.dumps(dict(export_script_sha256=sha(__file__),files=reduced),ensure_ascii=False,indent=2)+'\n')
        z.writestr('public_retrieval_manifest.json',json.dumps(dict(files=rows,full_native_movies_retained_on_ROG=excluded),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(archive=public.name,bytes=public.stat().st_size,sha256=sha(public),files=len(rows),native_movies_retained=len(excluded))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();export(a.root)
