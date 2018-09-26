set -e
cd gmx
rm forces.xvg
rm energy.xvg
gmx grompp -f martini_md.mdp -c minimized.gro -p dppc.top -o run.tpr
gmx mdrun -deffnm run -rerun minimized.gro -nt 1
echo pot | gmx energy -f run > energy.txt
echo 0 | gmx traj -f run -s run -of forces
rm mdout.mdp
rm run.edr run.log run.tpr run.trr
cd ..
