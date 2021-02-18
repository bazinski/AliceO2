<!-- doxy
\page refTRD TRD
/doxy -->


# DPL workflows for TRD

## Workflow and executables

* `o2-trd-trapsim` : trap chip simulation

* `o2-trd-tracklet-transformer` : calibrate the tracklets 

* `o2-trd-global-tracking` : build tracks @

* `o2-trd-trap2raw` : convert digit and/or tracklet data to raw data.

* `o2-trd-crucompressor` : Reads data in on the flp and re formulates it and sends it on to the stfsender.

* `o2-trd-epnreader` : Reads data in on the epn and formats it into messages and send it down the pipeline



```bash
# To decode digits from the raw (simulated) STF and send feed to to the workflow for further clusterization and reconstruction:
o2-raw-file-reader-workflow --loop 5 --delay 3 --conf ITSraw.cfg | o2-itsmft-stf-decoder-workflow --digits --no-clusters | o2-its-reco-workflow --disable-mc --digits-from-upstream
```

```bash
# To decode/clusterize the STF and feed directly the clusters into the workflow:
o2-raw-file-reader-workflow --loop 5 --delay 3 --conf ITSraw.cfg | o2-itsmft-stf-decoder-workflow | o2-its-reco-workflow --disable-mc --clusters-from-upstream
```

```bash
#If needed, one can request both digits and clusters from the STF decoder:
o2-raw-file-reader-workflow --loop 5 --delay 3 --conf ITSraw.cfg | o2-itsmft-stf-decoder-workflow --digits  | o2-its-reco-workflow --disable-mc --digits-from-upstream --clusters-from-upstream
```



## Tracking



## Associated devices :

* Digit Reader [o2::trd::TRDDigitReaderSpec](include/TRDWorkflow/TRDDigitReaderSpec.h)
* Trap Raw Writer [o2::trd::TRDTrapRawWriterSpec](include/TRDWorkflow/TRDTrapRawWriterSpec.h)
* Tracklet Writer [o2::trd::TRDTrackletWriterSpec](include/TRDWorkflow/TRDTrackletWriterSpec.h)
* Tracklet Reader [o2::trd::TRDTrackletReaderSpec](include/TRDWorkflow/TRDTrackletReaderSpec.h)
* Digit Writer [o2::trd::TRDDigitWriterSpec](include/TRDWorkflow/TRDDigitWriterSpec.h)
* Trap Simulator [o2::trd::TRDTrapSimulatorSpec](include/TRDWorkflow/TRDTrapSimulatorSpec.h)
* Digitizer [o2::trd::TRDDigitizerSpec](include/TRDWorkflow/TRDDigitizerSpec.h)
* Calibrated Tracklet Writer [o2::trd::TRDCalibratedTrackletWriterSpec](include/TRDWorkflow/TRDCalibratedTrackletWriterSpec.h) 
* Global Tracking [o2::trd::TRDGlobalTrackingSpec](include/TRDWorkflow/TRDGlobalTrackingSpec.h ) 
* TRDTrackingWorkflow.h    
* TRDTrackletTransformerSpec.h  
* TRDTrackWriterSpec.h     

