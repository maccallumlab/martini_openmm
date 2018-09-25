gmx grompp -f martini_md.mdp -c minimized.gro -p dppc.top -o run.tpr
gmx mdrun -deffnm run -rerun minimized.gro
