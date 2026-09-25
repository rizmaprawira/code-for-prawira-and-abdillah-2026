#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUTS_DIR="$BASE_DIR/inputs"
SCRIPT_TEMPLATE_SRC="/Users/rizzie/LinearBaroclinicModel/experiment/ncep_r2_2x2_moist_experiment/R1_p1BS_p1FRC/scripts"
RUN_EXP_TEMPLATE="/Users/rizzie/LinearBaroclinicModel/experiment/ncep_r2_2x2_moist_experiment/R1_p1BS_p1FRC/run_experiment.sh"

RUN_IDS=(
  CTRL
  EXP_B
  EXP_S-123
  EXP_B_S-123
  EXP_B_S-1
  EXP_B_S-12
)

BS_TAGS=(
  p1
  p2
  p1
  p2
  p2
  p2
)

FRC_TAGS=(
  p1
  p1
  p2
  p2
  p2_s1
  p2_s12
)

DESCRIPTIONS=(
  "Reference experiment with P1 basic state and P1 SST forcing over all regions."
  "Isolates basic-state effect: P2 basic state while retaining P1 SST forcing."
  "Isolates total SST-forcing effect: P2 forcing over all three regions with P1 basic state."
  "Full P2 configuration: P2 basic state and P2 SST forcing over all three regions."
  "Regional contribution: P2 SST forcing in Region 1 under P2 basic state (Regions 2 & 3 retain P1 forcing)."
  "Cumulative regional contribution: P2 SST forcing in Regions 1 & 2 under P2 basic state (Region 3 retains P1 forcing)."
)

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

require_file() {
  [[ -f "$1" ]] || die "Required source file is missing: $1"
}

atm_bs_years() {
  case "$1" in
    p1) printf '1981-2006\n' ;;
    p2) printf '2007-2025\n' ;;
    *) die "Unknown period tag: $1" ;;
  esac
}

echo "=== Staging 6 Comprehensive Moist LBM Experiments ==="

for index in "${!RUN_IDS[@]}"; do
  run_id="${RUN_IDS[$index]}"
  bs_tag="${BS_TAGS[$index]}"
  frc_tag="${FRC_TAGS[$index]}"
  desc="${DESCRIPTIONS[$index]}"
  target="$BASE_DIR/$run_id"
  bs_years="$(atm_bs_years "$bs_tag")"

  echo "Setting up $run_id (BS: $bs_tag, FRC: $frc_tag)..."
  mkdir -p "$target/bs/gt3/$bs_tag" "$target/bs/grads/$bs_tag" "$target/frc/$frc_tag" \
           "$target/config" "$target/scripts" "$target/out" "$target/post" "$target/log" "$target/plots"

  # Source files from inputs/
  atm_gt3="$INPUTS_DIR/basic_state/$bs_tag/ncep_r2_djf_${bs_years}.t21l20"
  atm_grads="$INPUTS_DIR/basic_state/$bs_tag/ncep_r2_djf_${bs_years}.t21l20.grd"
  sst_gt3="$INPUTS_DIR/basic_state/$bs_tag/sstwin.t21"
  sst_grads="$INPUTS_DIR/basic_state/$bs_tag/sstwin.t21.grd"
  sst_ctl="$INPUTS_DIR/basic_state/$bs_tag/sstwin.t21.ctl"

  topog="$INPUTS_DIR/basic_state/shared/grz.t21"
  lndidx="$INPUTS_DIR/basic_state/shared/gridx.t21"
  wg_gt3="$INPUTS_DIR/basic_state/shared/wgwin.t21"
  wg_grads="$INPUTS_DIR/basic_state/shared/wgwin.t21.grd"
  wg_ctl="$INPUTS_DIR/basic_state/shared/wgwin.t21.ctl"

  # Forcing file
  forcing="$INPUTS_DIR/forcing/$frc_tag/frcsst.${frc_tag}.t21.grd"
  forcing_ctl="$INPUTS_DIR/forcing/$frc_tag/frcsst.${frc_tag}.t21.ctl"

  for sf in "$atm_gt3" "$atm_grads" "$sst_gt3" "$sst_grads" "$sst_ctl" \
            "$topog" "$lndidx" "$wg_gt3" "$wg_grads" "$wg_ctl" "$forcing"; do
    require_file "$sf"
  done

  forcing_size="$(wc -c < "$forcing" | tr -d '[:space:]')"
  [[ "$forcing_size" == "16400" ]] || die "Forcing file must be 16400 bytes, found $forcing_size at $forcing"

  # Install inputs into run directory
  install -m 0644 "$atm_gt3" "$target/bs/gt3/$bs_tag/atm_basic_state.t21l20"
  install -m 0644 "$atm_grads" "$target/bs/grads/$bs_tag/atm_basic_state.t21l20.grd"
  install -m 0644 "$sst_gt3" "$target/bs/gt3/$bs_tag/sstwin.t21"
  install -m 0644 "$sst_grads" "$target/bs/grads/$bs_tag/sstwin.t21.grd"
  install -m 0644 "$sst_ctl" "$target/bs/grads/$bs_tag/sstwin.t21.ctl"

  install -m 0644 "$topog" "$target/bs/gt3/grz.t21"
  install -m 0644 "$lndidx" "$target/bs/gt3/gridx.t21"
  install -m 0644 "$wg_gt3" "$target/bs/gt3/wgwin.t21"
  install -m 0644 "$wg_grads" "$target/bs/grads/wgwin.t21.grd"
  install -m 0644 "$wg_ctl" "$target/bs/grads/wgwin.t21.ctl"

  install -m 0644 "$forcing" "$target/frc/$frc_tag/frcsst.${frc_tag}.t21.grd"
  if [[ -f "$forcing_ctl" ]]; then
    install -m 0644 "$forcing_ctl" "$target/frc/$frc_tag/frcsst.${frc_tag}.t21.ctl"
  fi

  # Copy scripts template
  cp -r "$SCRIPT_TEMPLATE_SRC"/* "$target/scripts/"

  # Generate config/experiment.env
  cat << ENV_EOF > "$target/config/experiment.env"
RUN_ID=${run_id}
BS_TAG=${bs_tag}
FRC_TAG=${frc_tag}
DESCRIPTION="${desc}"
ATM_BS_LOCAL=bs/gt3/${bs_tag}/atm_basic_state.t21l20
ATM_BS_GRADS_LOCAL=bs/grads/${bs_tag}/atm_basic_state.t21l20.grd
SSTFRC_LOCAL=frc/${frc_tag}/frcsst.${frc_tag}.t21.grd
SSTBS_LOCAL=bs/gt3/${bs_tag}/sstwin.t21
WGBS_LOCAL=bs/gt3/wgwin.t21
TOPOG_LOCAL=bs/gt3/grz.t21
LNDIDX_LOCAL=bs/gt3/gridx.t21
RUN_DAYS=30
ENV_EOF

  # Generate run_moist_lbm_generated.csh
  cat << CSH_EOF > "$target/run_moist_lbm_generated.csh"
#!/bin/csh -f
# Moist execution runner for ${run_id}

set ROOT = "${target}"
setenv LNHOME "/Users/rizzie/LinearBaroclinicModel/ln_solver"
setenv SYSTEM "mac"
setenv RUN "/Users/rizzie/LinearBaroclinicModel/ln_solver/model/bin/mac/lbm2.t21ml20ctintgr"
setenv DIR "\$ROOT/out"
setenv RSTFILE "Restart.amat"
setenv DATZ "../bs/gt3/grz.t21"
setenv DATS "../bs/gt3/${bs_tag}/sstwin.t21"
setenv DATW "../bs/gt3/wgwin.t21"
setenv DATI "../bs/gt3/gridx.t21"
setenv BSFILE "../bs/gt3/${bs_tag}/atm_basic_state.t21l20"
setenv FRC "../frc/${frc_tag}/frcsst.${frc_tag}.t21.grd"
setenv SFRC "../frc/${frc_tag}/frcsst.${frc_tag}.t21.grd"

cd "\$DIR"
echo job started at \`date\` > SYSOUT
/bin/rm -f SYSIN

cat << END_OF_DATA >! SYSIN
 &nmrun  run='moist linear model: ${run_id}'                     &end
 &nmtime start=0,1,1,0,0,0, end=0,1,30,0,0,0           &end
 &nmhdif order=4, tefold=6, tunit='HOUR'                       &end
 &nmdelt delt=40, tunit='MIN', inistp=2                        &end
 &nmdamp
  ddragv=1,1,1,5,15,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,
  ddragd=1,1,1,5,15,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,
  ddragt=1,1,1,5,15,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,
  tunit='DAY'
 &end
 &nminit file='\$BSFILE', DTBFR=0., DTAFTR=0., TUNIT='DAY'     &end
 &nmrstr file='\$RSTFILE', tintv=1, tunit='MON', overwt=t      &end
 &nmvdif vdifv=1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,
         vdifd=1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,
         vdift=1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3 &end
 &nmbtdif tdmpc=0.                                             &end
 &nmfrc  ffrc='\$FRC', oper=f, nfcs=1                         &end
 &nmsfrc fsfrc='\$SFRC', ofrc=t, osstf=t, nsfcs=1,
         fsend=-1,1,30,0,0,0 &end
 &nmmca ocum=t, ttau0=3, qtau0=3, ttauc=6, tunit='HOUR', sigkb=0.9D0 &end
 &nmsfcm owes=t, expw=0.5                                     &end
 &nmlsc olsc=t, ttaul0=2, qtaul0=2, tunit='HOUR', dqrat=1.D-2 &end
 &nmdata item='GRZ',   file='\$DATZ'                           &end
 &nmdata item='GRIDX', file='\$DATI'                           &end
 &nmdata item='GRWG',  file='\$DATW'                           &end
 &nmdata item='GRSST', file='\$DATS'                           &end
 &nmchck ocheck=f, ockall=f                                    &end
 &nmhisd tintv=1, tavrg=1, tunit='DAY'                         &end
 &nmhist item='PSI',  file='psi', tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='CHI',  file='chi', tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='U',    file='u',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='V',    file='v',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='OMGF', file='w',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='T',    file='t',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='Z',    file='z',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='PS',   file='p',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='Q',    file='q',   tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='Q1C',  file='dtc', tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='Q2C',  file='dqc', tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='Q1L',  file='dtl', tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='Q2L',  file='dql', tintv=1, tavrg=1, tunit='DAY' &end
 &nmhist item='PR',   file='pr',  tintv=1, tavrg=1, tunit='DAY' &end
END_OF_DATA

"\$RUN" < SYSIN |& tee -a SYSOUT
set run_status = \$status
echo job ended at \`date\` >> SYSOUT
exit \$run_status
CSH_EOF
  chmod +x "$target/run_moist_lbm_generated.csh"

  # Install run_experiment.sh
  install -m 0755 "$RUN_EXP_TEMPLATE" "$target/run_experiment.sh"

  echo "  Done staging $run_id"
done

echo "=== All 6 Experiments Staged Successfully ==="
