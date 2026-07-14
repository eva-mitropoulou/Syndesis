set ray_opaque_background, off
set antialias, 2
set ray_shadows, 0
set stick_radius, 0.18
set cartoon_transparency, 0.35
bg_color white

load results/md/reproducibility/mdcand_002/complex.gro, stable_system
hide everything, stable_system
show cartoon, stable_system and polymer
color lightblue, stable_system and polymer
select stable_ligand, stable_system and resn UNL
select stable_pocket, byres (stable_system and polymer within 5 of stable_ligand)
show sticks, stable_pocket
color gray70, stable_pocket
show sticks, stable_ligand
color teal, stable_ligand
orient stable_system
zoom stable_system, 2
ray 1800, 1350
png figures/manuscript/md_pocket_stable.png, 1600, 1200, 300, 1
delete all

load results/md/reproducibility/mdcand_neg01/complex.gro, negative_system
hide everything, negative_system
show cartoon, negative_system and polymer
color lightblue, negative_system and polymer
select negative_ligand, negative_system and resn UNL
select negative_pocket, byres (negative_system and polymer within 5 of negative_ligand)
show sticks, negative_pocket
color gray70, negative_pocket
show sticks, negative_ligand
color red, negative_ligand
orient negative_system
zoom negative_system, 2
ray 1800, 1350
png figures/manuscript/md_pocket_negative.png, 1600, 1200, 300, 1
quit
