set ray_opaque_background, off
set antialias, 2
set ray_shadows, 0
set stick_radius, 0.22
set cartoon_transparency, 0.25
bg_color white

load figures/manuscript/md_trajectory_frames/control_002_rep01_start.pdb, stable_start
load figures/manuscript/md_trajectory_frames/control_002_rep01_end.pdb, stable_end
align stable_end and polymer, stable_start and polymer

hide everything
show cartoon, stable_start and polymer
color lightblue, stable_start and polymer
select stable_start_ligand, stable_start and resn UNL
select stable_start_pocket, byres (stable_start and polymer within 5 of stable_start_ligand)
show sticks, stable_start_pocket
color gray70, stable_start_pocket
show sticks, stable_start_ligand
color teal, stable_start_ligand
orient stable_start and polymer
zoom stable_start and polymer, 1.03
ray 2200, 1350
png figures/manuscript/md_control_002_rep01_start.png, 2200, 1350, 300, 1

hide everything
show cartoon, stable_end and polymer
color lightblue, stable_end and polymer
select stable_end_ligand, stable_end and resn UNL
show sticks, stable_start_pocket
color gray70, stable_start_pocket
show sticks, stable_start_ligand
color gray80, stable_start_ligand
set stick_transparency, 0.55, stable_start_ligand
show sticks, stable_end_ligand
color teal, stable_end_ligand
ray 2200, 1350
png figures/manuscript/md_control_002_rep01_end.png, 2200, 1350, 300, 1
delete all

load figures/manuscript/md_trajectory_frames/misdocked_control_rep01_start.pdb, failed_start
load figures/manuscript/md_trajectory_frames/misdocked_control_rep01_end.pdb, failed_end
align failed_end and polymer, failed_start and polymer

hide everything
show cartoon, failed_start and polymer
color lightblue, failed_start and polymer
select failed_start_ligand, failed_start and resn UNL
select failed_start_pocket, byres (failed_start and polymer within 5 of failed_start_ligand)
show sticks, failed_start_pocket
color gray70, failed_start_pocket
show sticks, failed_start_ligand
color red, failed_start_ligand
orient failed_start and polymer
zoom failed_start and polymer, 1.03
ray 2200, 1350
png figures/manuscript/md_misdocked_control_rep01_start.png, 2200, 1350, 300, 1

hide everything
show cartoon, failed_end and polymer
color lightblue, failed_end and polymer
select failed_end_ligand, failed_end and resn UNL
show sticks, failed_start_pocket
color gray70, failed_start_pocket
show sticks, failed_start_ligand
color gray80, failed_start_ligand
set stick_transparency, 0.55, failed_start_ligand
show sticks, failed_end_ligand
color red, failed_end_ligand
ray 2200, 1350
png figures/manuscript/md_misdocked_control_rep01_end.png, 2200, 1350, 300, 1
quit
