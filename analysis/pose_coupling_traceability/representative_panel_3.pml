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
png representative_panel_3.png, 1000, 800, dpi=220, ray=1
