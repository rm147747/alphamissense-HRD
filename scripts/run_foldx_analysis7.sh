#!/bin/bash
# ============================================================
# AlphaHRD Analysis 7: FoldX ΔΔG for 27 PPI Interface Variants
# ============================================================
# Requirements:
#   - FoldX 5.0+ binary (academic license: https://foldxsuite.crg.eu/)
#   - Boltz-2 PDB files in structures/boltz_output/
#   - Python 3.8+ with pandas
#
# Usage:
#   chmod +x run_foldx_analysis7.sh
#   ./run_foldx_analysis7.sh /path/to/foldx
# ============================================================

FOLDX_BIN=${1:-"foldx"}
OUTDIR="foldx_results"
PDBDIR="structures/boltz_output"

mkdir -p $OUTDIR

echo "============================================================"
echo " FoldX ΔΔG Analysis — 27 Interface Variants"
echo "============================================================"

# Check FoldX
if ! command -v $FOLDX_BIN &>/dev/null; then
    echo "ERROR: FoldX not found at $FOLDX_BIN"
    echo "  Download from: https://foldxsuite.crg.eu/"
    echo "  Usage: $0 /path/to/foldx"
    exit 1
fi

# Step 1: Repair PDBs (optimize side chains)
echo ""
echo "[STEP 1] RepairPDB..."
for pdb in $PDBDIR/**/boltz_results_*/**/*_model_0.pdb; do
    name=$(basename $pdb .pdb)
    echo "  Repairing: $name"
    $FOLDX_BIN --command=RepairPDB --pdb=$pdb --output-dir=$OUTDIR 2>/dev/null
done

# Step 2: Generate individual_list.txt files per complex
echo ""
echo "[STEP 2] Generating mutation lists..."

# BRCA1_BARD1 complex mutations
cat > $OUTDIR/mutations_BRCA1_BARD1.txt << 'EOF'
CA24Y,A;
PA37L,A;
CA42R,A;
CA61G,A;
RA64C,A;
CY33B,B;
VM53B,B;
EOF

# BRCA1_PALB2 complex mutations
cat > $OUTDIR/mutations_BRCA1_PALB2.txt << 'EOF'
LA1396P,A;
EA1408K,A;
LP14B,B;
KR18B,B;
EOF

# PALB2_BRCA2 complex mutations
cat > $OUTDIR/mutations_PALB2_BRCA2.txt << 'EOF'
LA857P,A;
WC866A,A;
AT881A,A;
EOF

# BRCA2_RAD51 complex mutations
cat > $OUTDIR/mutations_BRCA2_RAD51.txt << 'EOF'
TA1523A,A;
DN1532A,A;
AT1548A,A;
EOF

# ATM_NBN complex mutations
cat > $OUTDIR/mutations_ATM_NBN.txt << 'EOF'
RC2644A,A;
SF2656A,A;
RH2700A,A;
RC740B,B;
DN745B,B;
EOF

# MRE11_NBN complex mutations
cat > $OUTDIR/mutations_MRE11_NBN.txt << 'EOF'
HY123A,A;
RC156A,A;
DN178A,A;
LP642B,B;
RH648B,B;
EOF

# Step 3: Run BuildModel for each complex
echo ""
echo "[STEP 3] Running FoldX BuildModel..."
for complex in BRCA1_BARD1 BRCA1_PALB2 PALB2_BRCA2 BRCA2_RAD51 ATM_NBN MRE11_NBN; do
    repaired_pdb="${OUTDIR}/${complex}_model_0_Repair.pdb"
    mutfile="${OUTDIR}/mutations_${complex}.txt"
    
    if [ -f "$repaired_pdb" ] && [ -f "$mutfile" ]; then
        echo "  Processing: $complex"
        $FOLDX_BIN --command=BuildModel \
            --pdb=$repaired_pdb \
            --mutant-file=$mutfile \
            --output-dir=$OUTDIR \
            --numberOfRuns=3 \
            2>/dev/null
    else
        echo "  SKIP: $complex (missing PDB or mutation file)"
        # Try with original PDB
        orig_pdb=$(find $PDBDIR -name "${complex}_model_0.pdb" 2>/dev/null | head -1)
        if [ -n "$orig_pdb" ]; then
            echo "    Trying original PDB: $orig_pdb"
            $FOLDX_BIN --command=BuildModel \
                --pdb=$orig_pdb \
                --mutant-file=$mutfile \
                --output-dir=$OUTDIR \
                --numberOfRuns=3 \
                2>/dev/null
        fi
    fi
done

# Step 4: Parse results
echo ""
echo "[STEP 4] Parsing ΔΔG results..."
python3 << 'PYEOF'
import os, glob, pandas as pd

outdir = "foldx_results"
results = []

for f in glob.glob(f"{outdir}/Dif_*.fxout") + glob.glob(f"{outdir}/Average_*.fxout"):
    try:
        df = pd.read_csv(f, sep='\t', skiprows=8)
        for _, row in df.iterrows():
            results.append({
                'file': os.path.basename(f),
                'mutation': row.get('Pdb', ''),
                'ddG': row.get('total energy', row.get('Total', None)),
            })
    except:
        pass

if results:
    df_res = pd.DataFrame(results)
    df_res.to_csv(f"{outdir}/foldx_ddg_summary.csv", index=False)
    
    print(f"{'Mutation':<30} {'ΔΔG (kcal/mol)':>15} {'Classification':>20}")
    print("-"*70)
    for _, r in df_res.iterrows():
        ddg = r['ddG']
        if ddg is not None:
            cls = "Neutral" if abs(ddg) < 1.0 else ("Destabilizing" if ddg > 1.0 else "Stabilizing")
            if abs(ddg) > 2.0: cls = "HIGHLY " + cls
            print(f"  {r['mutation']:<28} {ddg:>15.2f} {cls:>20}")
    
    n_destab = sum(1 for _, r in df_res.iterrows() if r['ddG'] and r['ddG'] > 1.0)
    print(f"\nTotal: {len(df_res)} variants scored")
    print(f"Destabilizing (ΔΔG > 1.0): {n_destab} ({100*n_destab/len(df_res):.0f}%)")
    print(f"Highly destabilizing (ΔΔG > 2.0): {sum(1 for _, r in df_res.iterrows() if r['ddG'] and r['ddG'] > 2.0)}")
else:
    print("No FoldX output files found. Check that FoldX ran successfully.")
    print("Expected output files: Dif_*.fxout or Average_*.fxout")
PYEOF

echo ""
echo "============================================================"
echo " DONE. Results in: $OUTDIR/foldx_ddg_summary.csv"
echo "============================================================"
