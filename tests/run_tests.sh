set -e

echo
echo
echo "Running simple_lipid test"
echo
echo
cd simple_lipid
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo
echo "Running complex_lipid test"
echo
echo
cd complex_lipid
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo
echo "Running protein test"
echo
echo
cd protein
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo
echo "Running membrane_protein test"
echo
echo
cd membrane_protein
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo
echo "Running martini3 aqp protein with en test"
echo
echo
cd 1j4n_en_m3
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo
echo "Running martini3 ubq protein with en test"
echo
echo
cd 1ubq_en_m3
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo
echo "Running martini3 bpti protein with en test"
echo
echo
cd 1k6u_en_m3
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo "All tests completed succesfully"
