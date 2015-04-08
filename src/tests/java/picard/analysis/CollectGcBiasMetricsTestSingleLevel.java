package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by kbergin on 4/8/15 to test GcBias on 'all reads' for comparison to a single level collector on a downsampled to 1% version of Solexa-272221.bam.
 */
public class CollectGcBiasMetricsTestSingleLevel extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectGcBiasMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "gc_bias_metrics_single_test.bam");
        final File Soutfile   = File.createTempFile("test", ".gc_bias_summary_metrics");
        final File Doutfile = File.createTempFile("test", ".gc_bias_detail_metrics");
        final File pdf   = File.createTempFile("test", ".pdf");
        final String referenceFile = "/Users/kbergin/picard/Homo_sapiens_assembly19.fasta";
        final String accLevel = "ALL_READS";
        final int windowSize = 100;
        final double minGenFraction = 1.0E-5;
        final boolean biSulfiteSeq = false;
        final boolean assumeSorted = false;
        Soutfile.deleteOnExit();
        Doutfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + Doutfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "SUMMARY_OUTPUT=" + Soutfile.getAbsolutePath(),
                "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                "WINDOW_SIZE=" + windowSize,
                "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                "METRIC_ACCUMULATION_LEVEL=" + accLevel,
                "ASSUME_SORTED=" + assumeSorted
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<GcBiasSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(Soutfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals("All Reads")) { //ALL_READS level
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 4089538);
                Assert.assertEquals(metrics.ALIGNED_READS, 8116117);
                Assert.assertEquals(metrics.AT_DROPOUT, 1.475997);
                Assert.assertEquals(metrics.GC_DROPOUT, 2.753563);
            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }
}

