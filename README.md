# ToT-analysis
Analysis of the ALICE 3 tracker PID using ToT

Requires O2 from https://github.com/hfribert/AliceO2, branch `ToT-pixcharge`
(stores the charge of every cluster pixel). Run everything from the O2 environment.

## Simulation

```bash
export ALICE3_SIM_FIELD=ON
export ALICE3_MAGFIELD_MACRO=<O2 source>/Detectors/Upgrades/ALICE3/macros/ALICE3Field.C

# TF3 = TOF (only as material here); add A3IP to -m to include the beam pipe
o2-sim-serial-run5 -n 50 -g boxgen -m TRK FT3 TF3 --seed 1 \
  --configKeyValues "Diamond.width[0]=0.01;Diamond.width[1]=0.01;Diamond.width[2]=5;BoxGun.pdg=211;BoxGun.number=20;BoxGun.eta[0]=-3;BoxGun.eta[1]=3;BoxGun.prange[0]=0.1;BoxGun.prange[1]=5;BoxGun.phirange[0]=0;BoxGun.phirange[1]=360;IOTOFBase.segmentedInnerTOF=true;IOTOFBase.segmentedOuterTOF=true;IOTOFBase.enableForwardTOF=false;IOTOFBase.enableBackwardTOF=false"
```

## Digitization

```bash
o2-steer-colcontexttool -i o2sim,200000,50:50 --use-existing-kine \
  --with-vertices kNoVertex --orbits 3 --orbitsPerTF 32

# ALMIRA response file; without --local-response-file the APTS response from CCDB is used
o2-sim-digitizer-workflow -b --run --onlyDet TRK,FT3 \
  --interactionRate 200000 \
  --incontext collisioncontext.root \
  --local-response-file <path>/ALICE3OTResponseData_v2.root \
  --configKeyValues "TRKDigitizerParam.noisePerPixel=1e-7;TRKDigitizerParam.chargeThreshold=91;FT3DigitizerParam.noisePerPixel=1e-7;FT3DigitizerParam.chargeThreshold=91;HBFUtils.nHBFPerTF=3"
```

## Clusterization

```bash
o2-alice3-trk-reco-workflow -b
```

## Cluster export

`ExportClusters.C` writes the TRK and FT3 clusters into a flat ROOT file, with the
seed pixel charge and the MC truth.

```bash
root -l -b -q '<repo>/ExportClusters.C+("<simulation directory>", "clusters.root")'
```

Tree `clusters` (one entry per cluster):

| Branch | Meaning |
|---|---|
| `det` | 0 = TRK, 1 = FT3 |
| `layer` | TRK 0-10 (VD 0-2, ML/OT 3-10), FT3 disk 0-11 |
| `chipID` | O2 sensor index |
| `row`, `col` | seed pixel on the sensor |
| `qSeed` | seed pixel charge [e-] |
| `x`, `y`, `z` | global cluster position [cm] |
| `event`, `trackID` | MC particle (-1 for noise) |
| `pdg` | PDG code |

Tree `particles`: `event`, `trackID`, `pdg`, `primary`, `mother`, vertex `vx`, `vy`, `vz` [cm], momentum `px`, `py`, `pz` [GeV/c].
