set -e

echo "Running simple_lipid test"
cd simple_lipid
pwd
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo "Running complex_lipid test"
cd complex_lipid
../Scripts/run_gmx.sh
../Scripts/run_openmm.sh
python ../Scripts/compare.py
cd ..
echo

echo
echo "All tests completed succesfully"
