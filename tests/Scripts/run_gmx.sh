set -e
cd gmx
gmx grompp -f martini_md.mdp -c minimized.gro -p system.top -n system.ndx -o run.tpr -maxwarn 1
gmx mdrun -deffnm run -rerun minimized.gro -nt 1
echo pot | gmx energy -f run > energy.txt
echo 0 | gmx traj -f run -s run -of forces
rm mdout.mdp
rm run.edr run.log run.tpr run.trr
cd ..
