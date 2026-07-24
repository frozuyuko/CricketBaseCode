ipsae = "~/CricketBaseCode/structureprediction/ipsae_ptm.py"
af3plddt = "~/CricketBaseCode/structureprediction/af3_plddt_ptm.py"
boltz2plddt = "~/CricketBaseCode/structureprediction/boltz2_plddt_ptm.py"

python "$ipsae" <path_to_pae_json_file> <path_to_mmcif_file> 10 10
python "$af3plddt" --cutoff 10 --chains A_B --out-prefix 
                        --out-file [OUT_FILE]
                        cif_path conf_json
python "$boltz2plddt" [--cutoff CUTOFF] [--chains CHAINS] [--print-chains]
                           [--out-prefix OUT_PREFIX] [--out-file OUT_FILE]
                           cif_path plddt_path

#Filter ipsae from highest ipsae max 
#Filter plddt from balanced chain interface plddt 

