# get .exe after building code
cp ../../1D_DG_solver_01.exe ./
# initialize the vertex points
python InitializeVertexPoints.py
# optimize vertex points
python DG_MDO.py