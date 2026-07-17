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
png representative_panel_2.png, 1000, 800, dpi=220, ray=1
