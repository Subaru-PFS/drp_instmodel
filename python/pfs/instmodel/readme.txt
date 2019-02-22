Sim2D dump from CPL

* 2D interpolation can be stripped out
* jegSpots stuff has been fluid
  + Caching:
    - FITS file with 43000 HDUs (!)
    - Convert to numpy x,y pickles per fiber for ease of reading; saved somewhere which may be "wrong"
  + 600 fibers in photons from simulator; normalise to the total flux of the brightest provided spot
* Operation:
  + Build traces one fiber at a time (5x oversampled)
  + Add to F32 image (non-oversampled) that is "just above" the detector
  + Conversion of this image to detector image is a separate module
  + Use real LAM bias, dark to convert to detector image
  + No flats because we don't have any real data
  + No NIR detector effects
  + Have a single realisation of the sky, but don't actually add it in
    - RHL sky swindle: set the variance to the photon noise from the sky plus a fraction of the sky levels
* External dependencies:
  + afw
  + cloomis instrument geometry (should be switched to afw)
* Biggest missing pieces:
  + pfiConfig (DAMD-32?)
  + Inverse flux calibration (throughput)
* Input spectra:
  + Line list (for arcs, combs)
  + Given spectrum (flats, science, sky)
  + No application of throughput
    - photons = (F/uJy) 15.09/R m^-2 s^-1
  + Sean Johnson is generating some simulated objects
    - Model spectra with proper wavelength sampling (needs to be oversampled)
    - CAN'T just use SDSS spectra because you get the LSF twice
  + Stars to come from Tanaka-san/NAOJ
* Tickets to consider:
  + DAMD-24: what is a pfsSim object?
    - pfsSim object is the spectrum specification
    - log-Wavelength WCS implemented, but should be tossed?
  + INFRA-208: where do we put them?
* Writes PFFA*.fits
  + The first "F" stands for "fake"
  + The headers are crazy (pyfits' fault)
* To do:
  + SIM2D-88 (Compare ETC and SIM2D+DRP2D outputs)
    - Include throughputs
    - Use same input spectra
    - Bug fixes
  + Survey simulation
    - pfiConfig
    - Sean Johnson input spectra (DAMD-24, INFRA-208)
  + b1 detectorMap
    - Craig will get this out of Jim

* To run:
  + lsst@87b210901008:~/pfs/drp_instmodel (master *$%=) $ makeSim -d r1 -F oneSky --exptime 900 -o oneSky.fits
  + For flats: `--imagetyp=flat --exptime=15`, for Ne arcs: `--imagetyp=arc,Ne --exptime=10`

Continuum:
E = n.h.nu = dnu.F_nu.T.A
n.h = F_nu.T.A.dnu/nu
h = 6.63e-34 J.s
1 Jy = 10^-26 W/m^2/Hz ==> 1nJy = 10^-35 W/m^2/Hz
n = 10^-35 W/m^2/Hz / 6.63e-34 J.s . (F_nu/nJy).T . (A/m^2) . (dnu/nu)
  = 0.015 . (F_nu/nJy).(T/s).(A/m^2).(dnu/nu)

1 Jy = 10^-26 W/m^2/Hz = 10^-23 erg/cm^2/s/Hz

E = n.h.nu = F.T.A
n = F.T.A.lambda/h/c
  = (F / W/m^2) . (T / sec) . (A / m^2) . 10^-9 . (lambda/nm) / 6.63e-34 J.s / 3e8 m/s


Original spectrum is 5000 nJy @ 680 nm

675-685 nm --> indices 582:697
Original spectrum is 5000 nJy @ 680 nm
In that region, we get about 104 +/- 14 counts in 900 sec.
Dispersion is 0.086 nm/pixel --> R = 7900
Transmission is 0.29

Expect: photons = 0.015 . (F_nu/nJy).(T/s).(A/m^2)/R = 130
Got: 104, so we're missing about 25% of the flux.



lsst@824ef520e03f:~/pfs/SIMS $ echo lsst.obs.pfs.PfsMapper > 20181213/_mapper
lsst@824ef520e03f:~/pfs/SIMS $ ingestImages.py /home/lsst/pfs/SIMS/20181213 20181213-raw/PFFA*.fits

sqlite> select * from raw;
id|site|category|field|visit|ccd|filter|arm|spectrograph|dateObs|expTime|dataType|taiObs|pfsConfigId|slitOffset
1|F|A|BIAS|0|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
2|F|A|BIAS|1|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
3|F|A|BIAS|2|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
4|F|A|BIAS|3|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
5|F|A|BIAS|4|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
6|F|A|BIAS|5|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
7|F|A|BIAS|6|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
8|F|A|BIAS|7|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
9|F|A|BIAS|8|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
10|F|A|BIAS|9|1|r|r|1|2018-12-13|0.0|BIAS|2018-12-13|0|0.0
11|F|A|DARK|10|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
12|F|A|DARK|11|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
13|F|A|DARK|12|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
14|F|A|DARK|13|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
15|F|A|DARK|14|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
16|F|A|DARK|15|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
17|F|A|DARK|16|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
18|F|A|DARK|17|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
19|F|A|DARK|18|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
20|F|A|DARK|19|1|r|r|1|2018-12-13|900.0|DARK|2018-12-13|0|0.0
21|F|A|FLAT|20|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
22|F|A|FLAT|21|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
23|F|A|FLAT|22|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
24|F|A|FLAT|23|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
25|F|A|FLAT|24|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
26|F|A|FLAT|25|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
27|F|A|FLAT|26|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
28|F|A|FLAT|27|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
29|F|A|FLAT|28|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
30|F|A|FLAT|29|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
31|F|A|FLAT|30|1|r|r|1|2018-12-13|30.0|FLAT|2018-12-13|0|0.0
32|F|A|ARC|31|1|r|r|1|2018-12-13|1.0|ARC|2018-12-13|0|0.0
33|F|A|OBJECT|32|1|r|r|1|2018-12-13|900.0|OBJECT|2018-12-13|0|0.0

lsst@824ef520e03f:~/pfs/SIMS $ mkdir 20181213/CALIB/

lsst@824ef520e03f:~/pfs/SIMS $ constructBias.py 20181213 --calib 20181213/CALIB --rerun bias --cores 3 --job bias --id field=BIAS
lsst@824ef520e03f:~/pfs/SIMS $ ingestCalibs.py 20181213 --calib 20181213/CALIB --validity 3000 20181213/rerun/bias/BIAS/pfsBias-2018-12-13-0-r1.fits

lsst@824ef520e03f:~/pfs/SIMS $ constructDark.py 20181213 --calib 20181213/CALIB --rerun dark --cores 3 --job dark --id field=DARK
lsst@824ef520e03f:~/pfs/SIMS $ ingestCalibs.py 20181213 --calib 20181213/CALIB --validity 3000 20181213/rerun/dark/DARK/pfsDark-2018-12-13-0-r1.fits

lsst@824ef520e03f:~/pfs/SIMS $ ingestCalibs.py 20181213 --calib 20181213/CALIB --validity 3000 $OBS_PFS_DIR/pfs/camera/pfsDetectorMap-005833-r1.fits --mode=link

lsst@824ef520e03f:~/pfs/SIMS $ constructFiberFlat.py 20181213 --calib 20181213/CALIB --rerun flat --cores 3 --job flat --id field=FLAT
lsst@824ef520e03f:~/pfs/SIMS $ ingestCalibs.py 20181213 --calib 20181213/CALIB --validity 3000 20181213/rerun/flat/FLAT/pfsFiberFlat-2018-12-13-000020-r1.fits

lsst@824ef520e03f:~/pfs/SIMS $ constructFiberTrace.py 20181213 --calib 20181213/CALIB --rerun fiberTrace --cores 1 --job fiberTrace --id visit=21
lsst@824ef520e03f:~/pfs/SIMS $ ingestCalibs.py 20181213 --calib 20181213/CALIB --validity 3000 20181213/rerun/fiberTrace/FIBERTRACE/pfsFiberTrace-2018-12-13-000021-r1.fits

lsst@824ef520e03f:~/pfs/SIMS $ reduceArc.py 20181213 --calib=20181213/CALIB/ --rerun arc --id field=ARC
lsst@824ef520e03f:~/pfs/SIMS $ reduceExposure.py 20181213 --calib 20181213/CALIB --rerun exp --id field=OBJECT


price@MacBook:~/pfs/spt_ExposureTimeCalculator (tickets/SIM2D-88=) $ python scripts/run_etc.py @scripts/run_etc.defaults --MAG_FILE=22.0 --MR_MODE=no --EXP_TIME=900 --OUTFILE_SNC=out/etc-22mag-900sec.dat --OUTFILE_SNL=-
price@MacBook:~/pfs/spt_ExposureTimeCalculator (tickets/SIM2D-88=) $ python scripts/gen_sim_spec.py @scripts/gen_sim_spec.defaults --outDir=out --etcFile=out/etc-22mag-900sec.dat --MAG_FILE=22.0

lsst@824ef520e03f:~/pfs/drp_instmodel (tickets/SIM2D-88 *%=) $ python armStats.py ~/pfs/spt_ExposureTimeCalculator/out/pfsArm-000001-r1.fits 
Signal: 0.375579683732
Noise: 0.061388810061
S/N: 6.11804795301

lsst@824ef520e03f:~/pfs/drp_instmodel (tickets/SIM2D-88 *%=) $ python armStats.py ~/pfs/SIMS/20181213/rerun/exp/pfsArm/2018-12-17/v0000033/pfsArm-000033-r1.fits 
Signal: 78.6044
Noise: 14.5264184933
S/N: 5.41113635056
