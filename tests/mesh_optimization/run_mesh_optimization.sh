# get .exe after building code
cp ../../DG_solver.exe ./
# initialize the vertex points
python InitializeVertexPoints.py
# optimize vertex points
python DG_MDO.py