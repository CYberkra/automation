from pathlib import Path
import gprMax
scene=gprMax.Scene()
scene.add(gprMax.Title(name=Path(__file__).stem))
scene.add(gprMax.Domain(p1=(6,6,6)))
scene.add(gprMax.Discretisation(p1=(.05,.05,.05)))
scene.add(gprMax.TimeWindow(time=200e-9))
scene.add(gprMax.OMPThreads(n=8))
scene.add(gprMax.PMLThickness(thickness=10))
scene.add(gprMax.Waveform(wave_type="impulse",amp=1,freq=1,id="impulse"))
scene.add(gprMax.HertzianDipole(polarisation="x",p1=(3,2.35,3),waveform_id="impulse"))
scene.add(gprMax.Rx(p1=(3,3.65,3),id="measurement",outputs=["Ex"]))
scene.add(gprMax.Rx(p1=(3,3,3),id="inside",outputs=["Ex"]))
gprMax.run(scenes=[scene],n=1,gpu=[0],gpu_precision="double",subgrid=False,autotranslate=True,outputfile=Path(__file__).with_suffix(""))
