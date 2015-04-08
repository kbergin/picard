/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;
import picard.analysis.directed.GcBiasMetricsCollector;
import picard.metrics.GcBiasMetrics;

import java.io.File;
import java.text.NumberFormat;
import java.util.List;
import java.util.Set;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Tool to collect information about GC bias in the reads in a given BAM file. Computes" +
                " the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%" +
                " and counts the number of read starts in each GC bin.  What is output and plotted is" +
                " the \"normalized coverage\" in each bin - i.e. the number of reads per window normalized" +
                " to the average number of reads per window across the whole genome..\n",
        usageShort = "Collects information about GC bias in the reads in the provided SAM or BAM",
        programGroup = Metrics.class
)
public class CollectGcBiasMetrics extends SinglePassSamProgram {
    /** The location of the R script to do the plotting. */
    private static final File R_SCRIPT = new File("/Users/kbergin/picard","gcBias.R");

    // Usage and parameters

    @Option(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Option(shortName = "S", doc = "The text file to write summary metrics to.", optional = true)
    public File SUMMARY_OUTPUT;

    @Option(doc = "The size of windows on the genome that are used to bin reads.")
    public int WINDOW_SIZE = 100;

    @Option(doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Option(shortName = "BS", doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Option(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates GcBiasMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private GcBiasMetricsCollector multiCollector;
    private String saveHeader;

    /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        if (SUMMARY_OUTPUT != null) IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        saveHeader = header.getReadGroups().get(0).getLibrary();
        //Delegate actual collection to GcBiasMetricCollector
        multiCollector = new GcBiasMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), WINDOW_SIZE, IS_BISULFITE_SEQUENCED);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Loop over the reference and the reads and calculate the basic metrics
    ////////////////////////////////////////////////////////////////////////////
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        multiCollector.acceptRecord(rec, ref);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Synthesize the normalized coverage metrics and write it all out to a file
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void finish() {
        multiCollector.finish();
        final MetricsFile<GcBiasMetrics, Integer> file = getMetricsFile();
        MetricsFile<GcBiasDetailMetrics, ?> detailMetricsFile = getMetricsFile();
        final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);
        final List<GcBiasMetrics> gcBiasMetricsList = file.getMetrics();
        for(final GcBiasMetrics gcbm : gcBiasMetricsList){
            List<GcBiasDetailMetrics> gcDetailList = gcbm.DETAILS.getMetrics();
            for(final GcBiasDetailMetrics d : gcDetailList) {
                detailMetricsFile.addMetric(d);
            }
            summaryMetricsFile.addMetric(gcbm.SUMMARY);
        }
        detailMetricsFile.write(OUTPUT);
        summaryMetricsFile.write(SUMMARY_OUTPUT);

        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        //I will need to edit the R script to output a chart for each details set
        RExecutor.executeFromFile(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                SUMMARY_OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                String.valueOf(WINDOW_SIZE));
    }
}


