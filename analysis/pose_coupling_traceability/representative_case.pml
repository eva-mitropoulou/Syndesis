reinitialize
bg_color white
set ray_opaque_background, off
set cartoon_transparency, 0.15
# Panel A: maximum CNNscore
load representative_poses/receptors/4hjo_a_aq4_1001.pdbqt, panel_1_receptor
load representative_poses/4hjo_a_aq4_1001_437219.pdbqt, panel_1_ligand
hide everything, panel_1_receptor
show cartoon, panel_1_receptor
color gray80, panel_1_receptor
show sticks, panel_1_ligand
util.cbag panel_1_ligand
select panel_1_pocket, byres (panel_1_receptor within 5 of panel_1_ligand)
show sticks, panel_1_pocket
color slate, panel_1_pocket
zoom panel_1_ligand, 12
save representative_1_4hjo_a_aq4_1001.pse
# Panel B: maximum recall
load representative_poses/receptors/1xkk_a_fmm_91.pdbqt, panel_2_receptor
load representative_poses/1xkk_a_fmm_91_437219.pdbqt, panel_2_ligand
hide everything, panel_2_receptor
show cartoon, panel_2_receptor
color gray80, panel_2_receptor
show sticks, panel_2_ligand
util.cbag panel_2_ligand
select panel_2_pocket, byres (panel_2_receptor within 5 of panel_2_ligand)
show sticks, panel_2_pocket
color slate, panel_2_pocket
zoom panel_2_ligand, 12
save representative_2_1xkk_a_fmm_91.pse
# Panel C: coupled-score selection
load representative_poses/receptors/1xkk_a_fmm_91.pdbqt, panel_3_receptor
load representative_poses/1xkk_a_fmm_91_437219.pdbqt, panel_3_ligand
hide everything, panel_3_receptor
show cartoon, panel_3_receptor
color gray80, panel_3_receptor
show sticks, panel_3_ligand
util.cbag panel_3_ligand
select panel_3_pocket, byres (panel_3_receptor within 5 of panel_3_ligand)
show sticks, panel_3_pocket
color slate, panel_3_pocket
zoom panel_3_ligand, 12
save representative_3_1xkk_a_fmm_91.pse
