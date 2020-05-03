#
# write_mzml <- function(){
#   system.time(
#     {
#       file <- tempfile()
#       header <- data.frame(
#         seqNum = 1:length(s2),
#         acquisitionNum = 1:length(s2),
#         msLevel = rep(2, length(s2)),
#         polarity = rep(1, length(2)),
#         peaksCount = vapply(map(s2, as, 'data.frame'), nrow, 1L, USE.NAMES = F),
#         totIonCurrent = as.vector(sapply(map(s2, intensity), sum)),
#         retentionTime = rep(1,length(s2)),
#         basePeakMZ = rep(1,length(s2)),
#         basePeakIntensity = rep(1,length(s2)),
#         collisionEnergy = rep(1,length(s2)),
#         ionisationEnergy = rep(1,length(s2)),
#         lowMZ = rep(1,length(s2)),
#         highMZ = rep(1,length(s2)),
#         precursorScanNum = rep(1,length(s2)),
#         precursorMZ = rep(1,length(s2)),
#         precursorCharge = rep(1,length(s2)),
#         precursorIntensity = rep(1,length(s2)),
#         mergedScan = rep(1,length(s2)),
#         mergedResultScanNum = rep(1,length(s2)),
#         mergedResultStartScanNum = rep(1,length(s2)),
#         mergedResultEndScanNum = rep(1,length(s2)),
#         filterString = as.character(names(s2)),
#         spectrumId = as.character(names(s2)),
#         centroided = rep(T, length(s2)),
#         injectionTime = rep(1,length(s2)),
#         ionMobilityDriftTime = rep(1,length(s2)), stringsAsFactors = F
#       )
#
#       mzR::writeMSData(object = map(map(s2, as, 'data.frame'), as.matrix), file = file, outformat = 'mzml', header = header)
#       msnexp <- MSnbase::readMSData(file = file)
#       sim_mat <- MSnbase::compareSpectra(msnexp)
#     }
#   )
#
# }
