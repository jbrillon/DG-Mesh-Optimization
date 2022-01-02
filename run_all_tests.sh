# h-convergence
cd ./tests/h_convergence/
./run_h_convergence.sh
python check_h_convergence.py
# # p-convergence
# cd ./tests/p_convergence/
# ./run_p_convergence.sh
# mesh optimization
cd ./tests/mesh_optimization/
./run_mesh_optimization.sh
python plotDGMDO.py