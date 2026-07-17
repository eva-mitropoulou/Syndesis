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
png representative_panel_1.png, 1000, 800, dpi=220, ray=1
