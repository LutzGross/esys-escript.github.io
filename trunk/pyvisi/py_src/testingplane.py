import vtk

file_name = "../test/python/data_data/interior_3D.xml"
xmlReader = vtk.vtkXMLUnstructuredGridReader()
xmlReader.SetFileName(file_name)

plane = vtk.vtkPlane()
planeCut = vtk.vtkCutter()
planeCut.SetInput(xmlReader.GetOutput())
planeCut.SetCutFunction(plane)

out = xmlReader.GetOutput()
center = out.GetCenter()
#origin = out.GetOrigin()

planeCut.Update()
plane.Modified()
#plane.SetOrigin(center[0], center[1], center[2]+ 0.001)
plane.SetOrigin(0.0,0.0,0.000000000000000000000000000001)

#--------------------------------------------
#plane.SetNormal(-0.01,0.0,1.0)
plane.SetNormal(0.0,0.0,1.0)
#------------------------------------------

planeCut.Update()
plane.Modified()

cutMapper = vtk.vtkDataSetMapper()
cutMapper.SetInput(planeCut.GetOutput())

cutActor = vtk.vtkActor()
cutActor.SetMapper(cutMapper)

outline=vtk.vtkOutlineFilter()
outline.SetInput(xmlReader.GetOutput())

outlineMapper = vtk.vtkDataSetMapper()
outlineMapper.SetInput(outline.GetOutput())

outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0,0,0)

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.AddActor(cutActor)
ren.AddActor(outlineActor)

ren.SetBackground(1,1,1)
renWin.SetSize(400,400)

iren.Initialize()
renWin.Render()
iren.Start()

