#!/bin/csh -f
# Moist execution runner for EXP_B_S-1

set ROOT = "/Users/rizzie/LinearBaroclinicModel/experiment/ncep_r2_moist_experiment/EXP_B_S-1"
setenv LNHOME "/Users/rizzie/LinearBaroclinicModel/ln_solver"
setenv SYSTEM "mac"
setenv RUN "/Users/rizzie/LinearBaroclinicModel/ln_solver/model/bin/mac/lbm2.t21ml20ctintgr"
setenv DIR "$ROOT/out"
setenv RSTFILE "Restart.amat"
setenv DATZ "../bs/gt3/grz.t21"
setenv DATS "../bs/gt3/p2/sstwin.t21"
setenv DATW "../bs/gt3/wgwin.t21"
setenv DATI "../bs/gt3/gridx.t21"
setenv BSFILE "../bs/gt3/p2/atm_basic_state.t21l20"
setenv FRC "../frc/p2_s1/frcsst.p2_s1.t21.grd"
setenv SFRC "../frc/p2_s1/frcsst.p2_s1.t21.grd"

cd "$DIR"
echo job started at `date` > SYSOUT
/bin/rm -f SYSIN

cat << END_OF_DATA >! SYSIN
 &nmrun  run='moist linear model: EXP_B_S-1'                     &end
 &nmtime start=0,1,1,0,0,0, end=0,1,30,0,0,0           &end
 &nmhdif order=4, tefold=6, tunit='HOUR'                       &end
 &nmdelt delt=40, tunit='MIN', inistp=2                        &end
 &nmdamp
  ddragv=1,1,1,5,15,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,
  ddragd=1,1,1,5,15,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,
  ddragt=1,1,1,5,15,30,30,30,30,30,30,30,30,30,30,30,30,30,1,1,
  tunit='DAY'
 &end
 &nminit file='$BSFILE', DTBFR=0., DTAFTR=0., TUNIT='DAY'     &end
 &nmrstr file='$RSTFILE', tintv=1, tunit='MON', overwt=t      &end
 &nmvdif vdifv=1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,
         vdifd=1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,
         vdift=1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3,1.d3 &end
 &nmbtdif tdmpc=0.                                             &end
 &nmfrc  ffrc='$FRC', oper=f, nfcs=1                         &end
 &nmsfrc fsfrc='$SFRC', ofrc=t, osstf=t, nsfcs=1,
         fsend=-1,1,30,0,0,0 &end
 &nmmca ocum=t, ttau0=3, qtau0=3, ttauc=6, tunit='HOUR', sigkb=0.9D0 &end
 &nmsfcm owes=t, expw=0.5                                     &end
 &nmlsc olsc=t, ttaul0=2, qtaul0=2, tunit='HOUR', dqrat=1.D-2 &end
 &nmdata item='GRZ',   file='$DATZ'                           &end
 &nmdata item='GRIDX', file='$DATI'                           &end
 &nmdata item='GRWG',  file='$DATW'                           &end
 &nmdata item='GRSST', file='$DATS'                           &end
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

"$RUN" < SYSIN |& tee -a SYSOUT
set run_status = $status
echo job ended at `date` >> SYSOUT
exit $run_status
