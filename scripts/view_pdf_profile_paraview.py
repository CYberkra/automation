"""Run with pvpython to save a physical-scale view of a prepared PDF profile."""
import argparse
from pathlib import Path
from paraview.simple import (
    XMLImageDataReader, GetActiveViewOrCreate, Show, ColorBy,
    GetColorTransferFunction, GetScalarBar, Text, Render, SaveState, SaveScreenshot,
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
out = args.out.resolve()
state = out/'line9_pdf_geometry.pvsm'
screenshot = out/'paraview_pdf_preview.png'
if state.exists() or screenshot.exists():
    raise ValueError('Refuse to overwrite a ParaView preview')
reader = XMLImageDataReader(registrationName='Line9 PDF: distance and elevation',
                           FileName=[str(out/'line9_pdf_geometry.vti')])
reader.UpdatePipeline()
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1600, 650]
view.Background = [.97, .97, .97]
view.UseColorPaletteForBackground = 0
view.InteractionMode = '2D'
display = Show(reader, view)
display.Representation = 'Surface'
display.Ambient = 1.
display.Diffuse = 0.
ColorBy(display, ('CELLS', 'MaterialID'))
lut = GetColorTransferFunction('MaterialID')
lut.InterpretValuesAsCategories = 1
lut.Annotations = ['0', 'air', '1', 'silty clay', '2', 'mudstone', '3', 'sandstone']
lut.IndexedColors = [217/255, 217/255, 217/255, 244/255, 177/255, 131/255,
                     191/255, 144/255, 0., 127/255, 96/255, 0.]
lut.IndexedOpacities = [1., 1., 1., 1.]
display.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = 'Material'
bar.ComponentTitle = ''
bar.TitleColor = [.15, .15, .15]
bar.LabelColor = [.15, .15, .15]
bar.WindowLocation = 'Upper Right Corner'
bar.ScalarBarLength = .3
display.DataAxesGrid.GridAxesVisibility = 1
display.DataAxesGrid.XTitle = 'Distance (m)'
display.DataAxesGrid.YTitle = 'Elevation (m)'
display.DataAxesGrid.ZTitle = 'Thin direction (m)'
for name in ('XTitleColor', 'YTitleColor', 'ZTitleColor', 'XLabelColor', 'YLabelColor', 'ZLabelColor'):
    setattr(display.DataAxesGrid, name, [.15, .15, .15])
view.CameraPosition = [175., 440., 500.]
view.CameraFocalPoint = [175., 440., .05]
view.CameraViewUp = [0., 1., 0.]
view.CameraParallelProjection = 1
view.CameraParallelScale = 78.
view.OrientationAxesVisibility = 1
note = Text(registrationName='Geometry provenance')
note.Text = ('Line9 PDF contacts | distance 0-350 m | elevation crop 415-465 m\n'
             'Equal horizontal/vertical scale | 5 cm viewing grid | no simulation\n'
             'Mudstone electrical properties pending; other materials are research assumptions')
text_display = Show(note, view)
text_display.WindowLocation = 'Lower Left Corner'
text_display.FontSize = 14
text_display.Color = [.15, .15, .15]
Render(view)
SaveState(str(state))
SaveScreenshot(str(screenshot), view, ImageResolution=[1600, 650])
print('ParaView state and physical-scale screenshot saved', flush=True)
