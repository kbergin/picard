package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.analysis.GcBiasDetailMetrics;
import picard.analysis.GcBiasSummaryMetrics;
import picard.metrics.GcBiasMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.metrics.MultiLevelCollector;
import picard.metrics.PerUnitMetricCollector;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashMap;

/**
 * Created by kbergin on 3/23/15.
 */
public class GcBiasMetricsCollector extends MultiLevelCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> {
    // Histograms to track the number of windows at each GC, and the number of read starts
    // at windows of each GC. Need 101 to get from 0-100.
    public static int LASTCONTIG = SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
    public static final int WINDOWS = 101;
    public byte[] refBases;
    public final int windowSize;
    public final boolean bisulfite;
    final int[] windowsByGc = new int[WINDOWS];
    byte[] gc;

    public GcBiasMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords, final int windowSize, final boolean bisulfite) {
        this.windowSize = windowSize;
        this.bisulfite = bisulfite;
        setup(accumulationLevels, samRgRecords);
    }

    // We will pass gc[] with the DefaultPerRecordCollectorArgs passed to the record collectors
    // This method is called once Per samRecord
    @Override
    protected GcBiasCollectorArgs makeArg(SAMRecord rec, ReferenceSequence ref) {
        if (ref!=null) {
            //only do the recalculation of gc if current ref is different from last ref
            if (ref.getContigIndex() != LASTCONTIG) {
                refBases = ref.getBases();
                StringUtil.toUpperCase(refBases);
                final int refLength = refBases.length;
                final int lastWindowStart = refLength - windowSize;
                gc = calculateAllGcs(refBases, windowsByGc, lastWindowStart);
                LASTCONTIG = ref.getContigIndex();
            }
        }
        return new GcBiasCollectorArgs(rec, gc);
    }

    /** Make a GcBiasCollector with the given arguments */
    @Override
    protected PerUnitMetricCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitGcBiasMetricsCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {super.acceptRecord(rec, ref);}

    /** A Collector for individual GcBiasMetrics for a given SAMPLE or SAMPLE/LIBRARY or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels) */
    public class PerUnitGcBiasMetricsCollector implements PerUnitMetricCollector<GcBiasMetrics, Integer, GcBiasCollectorArgs> {
        /** what do I do here - I think I need a map to store the different levels of collection results??*/
        Map<String, GcObject> GcData = new HashMap<String, GcObject>();
        String sample = null;
        String library = null;
        String readGroup = null;

        public PerUnitGcBiasMetricsCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            String prefix = null;
            if (this.readGroup != null) {
                prefix = this.readGroup;
                GcData.put(prefix, new GcObject());
            } else if (this.library != null) {
                prefix = this.library;
                GcData.put(prefix, new GcObject());
            } else if (this.sample != null) {
                prefix = this.sample;
                GcData.put(prefix, new GcObject());
            } else {
                prefix = "All_Reads.";
                GcData.put(prefix, new GcObject());
            }
        }

        public void acceptRecord(final GcBiasCollectorArgs args) {
            final SAMRecord rec = args.getRec();
            final byte[] gc = args.getGc();
            String type = null;
            if (this.readGroup != null) {
                type = this.readGroup;
                addRead(GcData.get(type), gc, rec);
            } else if (this.library != null) {
                type = this.library;
                addRead(GcData.get(type), gc, rec);
            } else if (this.sample != null) {
                type = this.sample;
                addRead(GcData.get(type), gc, rec);
            } else {
                type = "All_Reads.";
                addRead(GcData.get(type), gc, rec);
            }
        }

        public void finish() {}

        /** Sums the values in an int[]. */
        private double sum(final int[] values) {
            final int length = values.length;
            double total = 0;
            for (int i = 0; i < length; ++i) {
                total += values[i];
            }

            return total;
        }

        public void addMetricsToFile(final MetricsFile<GcBiasMetrics, Integer> file) {
            for (final Map.Entry<String, GcObject> entry : GcData.entrySet()) {
                final GcObject gcCur = entry.getValue();
                final int[] readsByGc = gcCur.readsByGc;
                final long[] errorsByGc = gcCur.errorsByGc;
                final long[] basesByGc = gcCur.basesByGc;
                final int totalClusters = gcCur.totalClusters;
                final int totalAlignedReads = gcCur.totalAlignedReads;

                /**maybe I need to put something here to finish before it's added to the file..
                 For each level in map of level -> gcObj do these calculations*/
                final GcBiasMetrics metrics = new GcBiasMetrics();
                final double totalWindows = sum(windowsByGc);
                final double totalReads = sum(readsByGc);
                final double meanReadsPerWindow = totalReads / totalWindows;

                for (int i = 0; i < windowsByGc.length; ++i) {
                    if (windowsByGc[i] == 0) continue;

                    final GcBiasDetailMetrics m = new GcBiasDetailMetrics();
                    m.GC = i;
                    m.WINDOWS = windowsByGc[i];
                    m.READ_STARTS = readsByGc[i];
                    if (errorsByGc[i] > 0) m.MEAN_BASE_QUALITY = QualityUtil.getPhredScoreFromObsAndErrors(basesByGc[i], errorsByGc[i]);
                    m.NORMALIZED_COVERAGE = (m.READ_STARTS / (double) m.WINDOWS) / meanReadsPerWindow;
                    m.ERROR_BAR_WIDTH = (Math.sqrt(m.READ_STARTS) / (double) m.WINDOWS) / meanReadsPerWindow;

                    metrics.DETAILS.add(m);
                }

                // Synthesize the high level metrics
                final GcBiasSummaryMetrics s = new GcBiasSummaryMetrics();
                s.WINDOW_SIZE = windowSize;
                s.TOTAL_CLUSTERS = totalClusters;
                s.ALIGNED_READS = totalAlignedReads;
                calculateDropoutMetrics(metrics.DETAILS, s);
                metrics.SUMMARY = s;
                file.addMetric(metrics);
            }
        }
    }

    /** Calculcate all the GC values for all windows. */
    private byte[] calculateAllGcs(final byte[] refBases, final int[] windowsByGc, final int lastWindowStart) {
        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];
        final CalculateGcState state = new CalculateGcState();
        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + windowSize;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
            if (windowGc != -1) windowsByGc[windowGc]++;
        }
        return gc;
    }

    /**
     * Calculates GC as a number from 0 to 100 in the specified window. If the window includes
     * more than five no-calls then -1 is returned.
     */
    private int calculateGc(final byte[] bases, final int startIndex, final int endIndex, final CalculateGcState state) {
        if (state.init) {
            state.init = false;
            state.gcCount = 0;
            state.nCount = 0;
            for (int i = startIndex; i < endIndex; ++i) {
                final byte base = bases[i];
                if (base == 'G' || base == 'C') ++state.gcCount;
                else if (base == 'N') ++state.nCount;
            }
        } else {
            final byte newBase = bases[endIndex - 1];
            if (newBase == 'G' || newBase == 'C') ++state.gcCount;
            else if (newBase == 'N') ++state.nCount;

            if (state.priorBase == 'G' || state.priorBase == 'C') --state.gcCount;
            else if (state.priorBase == 'N') --state.nCount;
        }
        state.priorBase = bases[startIndex];
        if (state.nCount > 4) return -1;
        else return (state.gcCount * 100) / (endIndex - startIndex);
    }

    /** Keeps track of current GC calculation state. */
    class CalculateGcState {
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
    /** Calculates the Illumina style AT and GC dropout numbers. */
    private void calculateDropoutMetrics(final Collection<GcBiasDetailMetrics> details,
                                         final GcBiasSummaryMetrics summary) {
        // First calculate the totals
        double totalReads = 0;
        double totalWindows = 0;

        for (final GcBiasDetailMetrics detail : details) {
            totalReads += detail.READ_STARTS;
            totalWindows += detail.WINDOWS;
        }

        double atDropout = 0;
        double gcDropout = 0;

        for (final GcBiasDetailMetrics detail : details) {
            final double relativeReads = detail.READ_STARTS / totalReads;
            final double relativeWindows = detail.WINDOWS / totalWindows;
            final double dropout = (relativeWindows - relativeReads) * 100;

            if (dropout > 0) {
                if (detail.GC <= 50) atDropout += dropout;
                if (detail.GC >= 50) gcDropout += dropout;
            }
        }

        summary.AT_DROPOUT = atDropout;
        summary.GC_DROPOUT = gcDropout;
    }

    /**Keeps track of each level of GcCalculation*/
    class GcObject{
        int totalClusters;
        int totalAlignedReads;
        int[] readsByGc;
        long[] basesByGc;
        long[] errorsByGc;

    }

    private void addRead(final GcObject gcObj, final byte[] gc, final SAMRecord rec) {
        if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++gcObj.totalClusters;
        if (!rec.getReadUnmappedFlag()) {
            final int pos = rec.getReadNegativeStrandFlag() ? rec.getAlignmentEnd() - windowSize : rec.getAlignmentStart();
            ++gcObj.totalAlignedReads;
            if (pos > 0) {
                final int windowGc = gc[pos];
                if (windowGc >= 0) {
                    ++gcObj.readsByGc[windowGc];
                    gcObj.basesByGc[windowGc] += rec.getReadLength();
                    gcObj.errorsByGc[windowGc] +=
                            SequenceUtil.countMismatches(rec, refBases, bisulfite) +
                                    SequenceUtil.countInsertedBases(rec) + SequenceUtil.countDeletedBases(rec);
                }
            }
        }
    }
}

// Arguments that need to be calculated once per SAMRecord that are then passed to each PerUnitMetricCollector
// for the given record
class GcBiasCollectorArgs {
    private final byte[] gc;
    private final SAMRecord rec;
    public byte[] getGc() {return gc;}
    public SAMRecord getRec() {return rec;}
    public GcBiasCollectorArgs(final SAMRecord rec, final byte[] gc) {
        this.gc = gc;
        this.rec = rec;
    }
}
