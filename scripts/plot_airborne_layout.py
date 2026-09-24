"""Render the actual M00 input geometry; no solver or field data."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    d={line.split(':',1)[0]:line.split(':',1)[1].split() for line in a.input.read_text().splitlines() if line.startswith('#') and ':' in line}
    domain=list(map(float,d['#domain']))
    cell=list(map(float,d['#dx_dy_dz']))
    pml=list(map(int,d['#pml_cells']))
    src=list(map(float,d['#hertzian_dipole'][1:4]))
    rx=list(map(float,d['#rx'][:3]))
    if '#material' in d or '#box' in d: raise ValueError('This drawing labels M00 all-air only')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    fig,axes=plt.subplots(1,2,figsize=(11,6),layout='constrained')
    for ax,vertical in zip(axes,(0,2)):
        w=domain[1];h=domain[vertical]
        ax.add_patch(Rectangle((0,0),w,h,facecolor='#eeeeee',edgecolor='black'))
        low_y=pml[1]*cell[1];high_y=pml[4]*cell[1]
        low_v=pml[vertical]*cell[vertical];high_v=pml[vertical+3]*cell[vertical]
        ax.add_patch(Rectangle((low_y,low_v),w-low_y-high_y,h-low_v-high_v,facecolor='white',edgecolor='#888888',linestyle=':'))
        ax.plot([src[1],rx[1]],[src[vertical],rx[vertical]],color='#333333')
        ax.scatter(src[1],src[vertical],s=70,color='#0072B2',label='Tx: ideal x-current element',zorder=5)
        ax.scatter(rx[1],rx[vertical],s=70,marker='s',color='#D55E00',label='Rx: Ex sample',zorder=5)
        ax.annotate('Tx',(src[1],src[vertical]),xytext=(-23,12),textcoords='offset points')
        ax.annotate('Rx',(rx[1],rx[vertical]),xytext=(8,12),textcoords='offset points')
        ax.set(xlim=(-.5,w+.5),ylim=(-.5,h+.5),xlabel='Y: cross-track [m]',ylabel=('X: flight direction [m]' if vertical==0 else 'Z: height [m]'))
        ax.set_aspect('equal');ax.set_xticks(range(0,25,4));ax.set_yticks(range(0,25,4))
        ax.text(1.5,1.7,'Shaded rim: 1 m HORIPML',fontsize=8,color='#555555')
    axes[0].annotate('Flight +X',xy=(5,18),xytext=(5,9),arrowprops=dict(arrowstyle='->',lw=2),ha='center')
    axes[0].annotate('1.3 m transverse baseline',xy=(11.95,12),xytext=(12,7),arrowprops=dict(arrowstyle='-'),ha='center',fontsize=9)
    axes[0].set_title('Plan view: X–Y')
    axes[1].axhline(4,ls='--',color='#777777')
    axes[1].text(1.5,4.6,'Z=4 m: future ground reference only',fontsize=8)
    axes[1].annotate('',xy=(19,19),xytext=(19,4),arrowprops=dict(arrowstyle='<->'))
    axes[1].text(19.7,11.5,'15 m',va='center',fontsize=9)
    axes[1].set_title('Cross-track section: Y–Z')
    axes[0].legend(loc='upper left',fontsize=8)
    fig.suptitle('M00: entire domain is air; 5 cm cells, 400 ns requested record\nX polarisation is a modelling assumption, not confirmed antenna orientation',fontsize=11)
    fig.savefig(a.output_dir/'layout.png',dpi=160)
    plt.close(fig)
    (a.output_dir/'provenance.json').write_text(json.dumps(dict(input=str(a.input),input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),solver_called=False),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
