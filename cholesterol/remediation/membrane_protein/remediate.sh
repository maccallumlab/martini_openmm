gmx grompp -f minimize.mdp -c before_vsites.gro -p system.top -o minimize -n system.ndx -maxwarn 3
gmx_d mdrun -s minimize.tpr -c minimized_not_whole.gro -nt 1 -v
echo 0 | gmx trjconv -f minimized_not_whole.gro -s minimize.tpr -o minimized.gro -pbc whole
rm \#*