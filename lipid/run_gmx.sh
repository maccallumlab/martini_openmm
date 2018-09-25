cd gmx
gmx grompp -f martini_md.mdp -c minimized.gro -p dppc.top -o run.tpr
gmx mdrun -deffnm run -rerun minimized.gro
echo pot | gmx energy -f run > energy.txt
echo 0 | gmx traj -f run -s run -of forces
cd ..
