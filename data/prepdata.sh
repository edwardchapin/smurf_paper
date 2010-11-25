#!/bin/tcsh
#+
#  Script for preparing data for plots. Currently set up to run on octopus,
#  but you can fiddle dataroot which is the path to a tree resembling
# /jcmtdata/raw/scuba2 at JAC
#

set dataroot = /scuba2
set msg = quiet

# the following outputs are used by intro.pro ---------------------------------

echo "*** Processing 20091214_00015 ***"

echo "    Concatenating..."

$SMURF_DIR/sc2concat "$dataroot/s4a/20091214/00015/s4a20091214_00015_000[23]" \
    s4a20091214_00015_con msg_filter=$msg

$SMURF_DIR/sc2concat "$dataroot/s8d/20091214/00015/s8d20091214_00015_000[23]" \
    s8d20091214_00015_con msg_filter=$msg

echo "    Cleaning..."

$SMURF_DIR/sc2clean s4a20091214_00015_con s4a20091214_00015_con_clean \
    config=^$STARLINK_DIR/share/smurf/dimmconfig.lis msg_filter=$msg

$SMURF_DIR/sc2clean s8d20091214_00015_con s8d20091214_00015_con_clean \
    config=^$STARLINK_DIR/share/smurf/dimmconfig.lis msg_filter=$msg

echo "    Convert to FITS..."

if( -e s4a20091214_00015_con_clean.fits ) then
  rm s4a20091214_00015_con_clean.fits
endif

$CONVERT_DIR/ndf2fits s4a20091214_00015_con_clean \
    s4a20091214_00015_con_clean.fits msg_filter=$msg

if( -e s8d20091214_00015_con_clean.fits ) then
  rm s8d20091214_00015_con_clean.fits
endif

$CONVERT_DIR/ndf2fits s8d20091214_00015_con_clean \
    s8d20091214_00015_con_clean.fits msg_filter=$msg

echo "    Getting JCMTState..."

$SMURF_DIR/jcmtstate2cat s4a20091214_00015_con > state_20091214_00015.tst
