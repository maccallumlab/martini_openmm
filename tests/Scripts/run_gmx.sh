set -e
cd gmx
# use double precision for everything
gmx_d grompp -f martini_md.mdp -c minimized.gro -p system.top -n system.ndx -o run.tpr -r minimized.gro -maxwarn 1 -v
gmx_d mdrun -deffnm run -rerun minimized.gro -nt 1 -v
echo pot | gmx_d energy -f run > energy.txt
echo 0 | gmx_d traj -f run -s run -of forces
rm mdout.mdp
rm run.edr run.log run.tpr run.trr
cd ..
