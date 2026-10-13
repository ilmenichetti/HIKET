#!/bin/bash
# =============================================================================
# refresh_all_local.sh -- rebuild EVERYTHING downstream of a calibration, on the Mac (2026-10-13)
#
# The sequence of the full refresh of 2026-10-13, in order, as one command. Run after
# the arms have landed on Roihu and been synced (step 0, by hand -- CLAUDE.md):
#   rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/runs/        ./Calibration_real_data_transient/runs/
#   rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/diagnostics/ ./Calibration_real_data_transient/diagnostics/
#   rsync -av roihu:/scratch/project_2019134/HIKET/Data/model_inputs/                           ./Data/model_inputs/
#
#   1  production: predictive, residuals, multimodel comparison (stages 2-4)
#   2  balanced-start arm: predictive + F16-F19, T_eqinit
#   3  analyses: transit times, initial state, headroom, input to hold the sink,
#      plot-scale checks, the published arm and the forward experiment, and the
#      decomposition -- fixed parameters, balanced-start parameters, and the RECALIBRATED
#      version when its arms exist (_eqinit_climheld; five terms if _eqinit_litheld and
#      _climheld exist too)
#   4  every manuscript figure builder (F4 first: F3 and T1 read its cache)
#   5  the residual chain (Yasso07+15 target) and its two manuscript copies
#   6  the NextGenC bundle (matrices, maps, thumbnails, delta, vulnerability) + SOC_maps sync
#   7  clean the setup by-products: input bundles / diagnostics keyed to RUN_IDs that have
#      no posterior (written by every script that sources a calibration script)
#
# ⚠ Every script that builds a starting state must use the run's OWN pre-run climate
#   (hiket_set_prerun_clim from the input bundle). F15 did not, and was wrong for a day.
#
# Usage (repo root):  bash Calibration_real_data_transient/refresh_all_local.sh [--from N] [--skip-report]
# Logs:   $LOGDIR (default /tmp/hiket_refresh_<timestamp>); one log per script, summary at the end.
# Afterwards, by hand: restate the numbers (draft, storyline, appendices; NEXT_SESSION.md
# before/after table), compile the documents, commit.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."
FROM=1; REPORT=1
while [ $# -gt 0 ]; do case "$1" in --from) FROM=$2; shift 2 ;; --skip-report) REPORT=0; shift ;; *) echo "unknown: $1"; exit 1 ;; esac; done
LOGDIR=${LOGDIR:-/tmp/hiket_refresh_$(date +%Y%m%d_%H%M%S)}; mkdir -p "$LOGDIR"
MARK="$LOGDIR/.start"; touch "$MARK"
FAILS=()
run() {   # run <step-label> <command...>
  local lab=$1; shift; local s=$(date +%s)
  "$@" > "$LOGDIR/${lab}.log" 2>&1; local rc=$?
  printf "%-44s exit %d  (%ds)\n" "$lab" $rc $(( $(date +%s) - s )) | tee -a "$LOGDIR/summary.txt"
  [ $rc -eq 0 ] || FAILS+=("$lab")
}
R() { Rscript --no-save "$@"; }
has_arm() { ls Calibration_real_data_transient/runs/*"$1"_posterior_*.rds >/dev/null 2>&1; }
echo "logs: $LOGDIR" | tee "$LOGDIR/summary.txt"

if [ $FROM -le 1 ]; then
  run 1_production_stages_2-4 R Calibration_real_data_transient/run_hiket_pipeline.R --skip-calibration
fi
if [ $FROM -le 2 ]; then
  run 2_eqinit_predictive_and_figures bash manuscript/figures/run_eqinit_figures.sh
fi
if [ $FROM -le 3 ]; then
  for f in intrinsic_mrt init_state_plausibility equilibrium_headroom input_to_hold_sink basal_area_sign \
           residuals_by_latitude yasso20_late_source_check woody_compartment_test \
           published_arm_dat f15_forward_experiment attribution_decomposition; do
    run 3_$f R doublechecks/$f.R
  done
  run 3_attribution_eqparams env HIKET_ATTR_PARAMS=eq Rscript --no-save doublechecks/attribution_decomposition.R
  if has_arm _eqinit_climheld; then
    run 3_attribution_recalibrated env HIKET_ATTR_MODE=recalibrated Rscript --no-save doublechecks/attribution_decomposition.R
  else echo "3_attribution_recalibrated: SKIPPED (no _eqinit_climheld arm)" | tee -a "$LOGDIR/summary.txt"; fi
fi
if [ $FROM -le 4 ]; then
  run 4_build_F4_initialization R manuscript/figures/build_F4_initialization.R
  for f in manuscript/figures/build_*.R; do
    b=$(basename "$f" .R)
    case "$b" in build_F4_initialization|build_F16_*|build_F17_*|build_F18_*|build_F19_*|build_T_eqinit) continue ;; esac
    run 4_$b R "$f"
  done
  run 4_S15S17_joint env HIKET_ATTR_JOINT=1 Rscript --no-save manuscript/figures/build_S15S17_attribution.R
  if [ -f doublechecks/attribution_recalibrated.rds ] && has_arm _eqinit_climheld; then
    run 4_S15S17_recal env HIKET_ATTR_RECAL=1 Rscript --no-save manuscript/figures/build_S15S17_attribution.R
    if has_arm _eqinit_litheld && has_arm _climheld; then
      run 4_S15S17_recal_joint env HIKET_ATTR_RECAL=1 HIKET_ATTR_JOINT=1 Rscript --no-save manuscript/figures/build_S15S17_attribution.R
    fi
  fi
fi
if [ $FROM -le 5 ]; then
  for f in residuals_rf_everything residuals_stand_and_soil residuals_time_since_cut residuals_rf_top_scatter \
           residuals_basal_area_by_year residuals_tally_predictors residuals_chemistry_predictors \
           residuals_events_over_time build_F11_alt_rf_heatmap harvest_dynamics rh_vs_stand_age understorey_litter; do
    [ -f doublechecks/$f.R ] && run 5_$f R doublechecks/$f.R
  done
  cp doublechecks/figures/residuals_stand_and_soil.png      manuscript/figures/F11_stand_and_soil.png
  cp doublechecks/figures/residuals_time_since_cut_summary.png manuscript/figures/F11b_time_since_cut.png
fi
if [ $FROM -le 6 ] && [ $REPORT -eq 1 ]; then
  for f in build_soc_matrices plot_soc_trajectories build_soc_maps make_thumbnails make_delta_map \
           build_vulnerability_map diagnostic_semivariogram; do
    run 6_$f bash -c "cd Reporting/NextgenC_report && Rscript --no-save $f.R"
  done
  (cd Reporting/NextgenC_report && cp *_SOC_mean.csv *_SOC_sd.csv Finland_SOC_matrices.ods SOC_maps/) \
    && echo "6_SOC_maps_sync                              done" | tee -a "$LOGDIR/summary.txt"
fi

# 7  setup by-products: bundles / diagnostics keyed to a RUN_ID with no posterior, created by this refresh
n=0
for f in $(find Data/model_inputs Calibration_real_data_transient/diagnostics -newer "$MARK" -type f \
             -name "*_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]*" 2>/dev/null); do
  id=$(basename "$f" | grep -oE "[0-9]{8}_[0-9]{6}" | head -1)
  ls Calibration_real_data_transient/runs/*_posterior_${id}.rds >/dev/null 2>&1 || { rm -f "$f"; n=$((n+1)); }
done
echo "7_cleanup: removed $n setup by-products" | tee -a "$LOGDIR/summary.txt"

echo; if [ ${#FAILS[@]} -eq 0 ]; then echo "ALL STEPS OK"; else echo "FAILED: ${FAILS[*]}  (logs in $LOGDIR)"; exit 1; fi
