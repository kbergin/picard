package picard.analysis;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.DbSnpBitSetUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static java.lang.Math.log10;

/**
 * Quantify substitution errors caused by mismatched base pairings during various
 * stages of sample / library prep. These are most commonly caused by oxidation
 * (e.g. the 8-oxo-G error mode), but could have other causes as well.
 *
 * We measure two distinct error types - artifacts that are introduced before
 * the addition of palindromic adapters ("pre adapt") and those that are
 * introduced after the addition of bait targets ("bait bias"). For each of
 * these, we provide summary metrics as well as detail metrics broken down
 * by reference context.
 *
 * @author mattsooknah
 *
 */
@CommandLineProgramProperties(
        usage = CollectDnaOxidationMetrics.USAGE,
        usageShort = CollectDnaOxidationMetrics.USAGE,
        programGroup = Metrics.class
)
public class CollectDnaOxidationMetrics extends CommandLineProgram {
    static final String USAGE = "Collect metrics relating to various kinds of DNA oxidation damage.";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input BAM file for analysis.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Base path of output files to write.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME,
            doc = "Reference sequence to which BAM is aligned.")
    public File REFERENCE_SEQUENCE;

    @Option(doc = "An optional list of intervals to restrict analysis to.",
            optional = true)
    public File INTERVALS;

    @Option(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.",
            optional = true)
    public File DB_SNP;

    @Option(shortName = "Q",
            doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName = "MQ",
            doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName = "MIN_INS",
            doc = "The minimum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName = "MAX_INS",
            doc = "The maximum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    // TODO rename this to CONTEXT_RADIUS? CONTEXT_SIZE could be mistakenly interpreted as including both sides
    @Option(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "The optional set of sequence contexts to restrict analysis to. Their reverse complements will also be analyzed, " +
            "even if not specified here. If not supplied all contexts are analyzed.")
    public Set<String> CONTEXTS = new HashSet<String>();

    @Option(doc = "For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public int STOP_AFTER = Integer.MAX_VALUE;

    private final Log log = Log.getInstance(CollectDnaOxidationMetrics.class);
    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";
    private static final byte[] BASES = {'A', 'C', 'G', 'T'};
    private static final double MIN_ERROR = 1e-10; // minimum error rate to report

    // Stock main method
    public static void main(final String[] args) {
        new CollectDnaOxidationMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final int size = 1 + 2 * CONTEXT_SIZE;
        final List<String> messages = new ArrayList<String>();

        for (final String ctx : CONTEXTS) {
            if (ctx.length() != size) {
                messages.add("Context " + ctx + " is not " + size + " long as implied by CONTEXT_SIZE=" + CONTEXT_SIZE);
            }
        }

        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected int doWork() {
        final File PRE_ADAPT_SUMMARY_OUT = new File(OUTPUT + ".pre_adapt_summary_metrics");
        final File PRE_ADAPT_DETAILS_OUT = new File(OUTPUT + ".pre_adapt_detail_metrics");
        final File BAIT_BIAS_SUMMARY_OUT = new File(OUTPUT + ".bait_bias_summary_metrics");
        final File BAIT_BIAS_DETAILS_OUT = new File(OUTPUT + ".bait_bias_detail_metrics");

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(PRE_ADAPT_SUMMARY_OUT);
        IOUtil.assertFileIsWritable(PRE_ADAPT_DETAILS_OUT);
        IOUtil.assertFileIsWritable(BAIT_BIAS_SUMMARY_OUT);
        IOUtil.assertFileIsWritable(BAIT_BIAS_DETAILS_OUT);

        if (INTERVALS != null) IOUtil.assertFileIsReadable(INTERVALS);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);

        final Set<String> samples = new HashSet<String>();
        final Set<String> libraries = new HashSet<String>();
        for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
            samples.add(nvl(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(nvl(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        // Setup the calculator
        final int contextLength = (2 * CONTEXT_SIZE) + 1;
        final Set<String> contexts = CONTEXTS.isEmpty() ? makeContextStrings(contextLength) : includeReverseComplements(CONTEXTS);
        final ArtifactCounter counts = new ArtifactCounter(libraries, contexts);

        // Load up dbSNP if available
        log.info("Loading dbSNP File: " + DB_SNP);
        final DbSnpBitSetUtil dbSnp;
        if (DB_SNP != null) dbSnp = new DbSnpBitSetUtil(DB_SNP, in.getFileHeader().getSequenceDictionary());
        else dbSnp = null;

        // Make an iterator that will filter out funny looking things
        final SamLocusIterator iterator;
        if (INTERVALS != null) {
            final IntervalList intervals = IntervalList.fromFile(INTERVALS);
            iterator = new SamLocusIterator(in, intervals.uniqued(), false);
        } else {
            iterator = new SamLocusIterator(in);
        }
        iterator.setEmitUncoveredLoci(false);
        iterator.setMappingQualityScoreCutoff(MINIMUM_MAPPING_QUALITY);
        iterator.setSamFilters(Arrays.asList(
                new NotPrimaryAlignmentFilter(),
                new DuplicateReadFilter(),
                new InsertSizeFilter(MINIMUM_INSERT_SIZE, MAXIMUM_INSERT_SIZE)
        ));

        log.info("Starting iteration.");
        long nextLogTime = 0;
        int sites = 0;

        // Iterate over reference loci
        for (final SamLocusIterator.LocusInfo info : iterator) {
            // Skip dbSNP sites
            final String chrom = info.getSequenceName();
            final int pos = info.getPosition();
            final int index = pos - 1;
            if (dbSnp != null && dbSnp.isDbSnpSite(chrom, pos)) continue;

            // Skip sites at the end of chromosomes
            final byte[] bases = refWalker.get(info.getSequenceIndex()).getBases();
            if (pos < contextLength || pos > bases.length - contextLength) continue;

            // TODO we should account for SNPs/indels in the surrounding context

            // Get the reference context string and perform counting
            final String context = StringUtil.bytesToString(bases, index - CONTEXT_SIZE, contextLength).toUpperCase();
            if (contexts.contains(context)) counts.countAlleles(info, context);

            // See if we need to stop
            if (++sites % 100 == 0) {
                final long now = System.currentTimeMillis();
                if (now > nextLogTime) {
                    log.info("Visited " + sites + " sites of interest. Last site: " + chrom + ":" + pos);
                    nextLogTime = now + 60000;
                }
            }
            if (sites >= STOP_AFTER) break;
        }

        // Finish up and write metrics
        final MetricsFile<PreAdaptSummaryMetrics, Integer> preAdaptSummaryMetricsFile = getMetricsFile();
        final MetricsFile<PreAdaptDetailMetrics, Integer> preAdaptDetailMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasSummaryMetrics, Integer> baitBiasSummaryMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasDetailMetrics, Integer> baitBiasDetailMetricsFile = getMetricsFile();

        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        final List<LibraryLevelMetrics> allMetrics = counts.finalize(sampleAlias);

        for (final LibraryLevelMetrics llm : allMetrics) {
            for (final PreAdaptDetailMetrics padm : llm.preAdaptDetailMetrics) preAdaptDetailMetricsFile.addMetric(padm);
            for (final BaitBiasDetailMetrics bbdm : llm.baitBiasDetailMetrics) baitBiasDetailMetricsFile.addMetric(bbdm);
            preAdaptSummaryMetricsFile.addMetric(llm.preAdaptSummaryMetrics);
            baitBiasSummaryMetricsFile.addMetric(llm.baitBiasSummaryMetrics);
        }

        preAdaptDetailMetricsFile.write(PRE_ADAPT_DETAILS_OUT);
        preAdaptSummaryMetricsFile.write(PRE_ADAPT_SUMMARY_OUT);
        baitBiasDetailMetricsFile.write(BAIT_BIAS_DETAILS_OUT);
        baitBiasSummaryMetricsFile.write(BAIT_BIAS_SUMMARY_OUT);

        CloserUtil.close(in);
        return 0;
    }

    // --------------------------------------------------------------------------------------------
    // Core extraction methods
    // --------------------------------------------------------------------------------------------

    /**
     * Core extraction method for detail metrics.
     */
    private DetailMetricsCollector extractDetailMetrics(final String sampleAlias, final String library, final ContextAccumulator contextAccumulator) {
        final DetailMetricsCollector detailMetrics = new DetailMetricsCollector();
        for (final String context : contextAccumulator.keySet()) {
            final byte refBase = (byte) getCentralBase(context);
            for (final byte altBase : BASES) {
                if (altBase != refBase) {
                    final PreAdaptDetailMetrics padm = new PreAdaptDetailMetrics();
                    final BaitBiasDetailMetrics bbdm = new BaitBiasDetailMetrics();

                    // retrieve all the necessary alignment counters
                    final AlignmentAccumulator fwdRefAlignments = contextAccumulator.get(context).get(refBase);
                    final AlignmentAccumulator fwdAltAlignments = contextAccumulator.get(context).get(altBase);
                    final AlignmentAccumulator revRefAlignments = contextAccumulator.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(refBase));
                    final AlignmentAccumulator revAltAlignments = contextAccumulator.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(altBase));

                    // populate basic fields
                    padm.SAMPLE_ALIAS = sampleAlias;
                    padm.LIBRARY = library;
                    padm.CONTEXT = context;
                    padm.REF_BASE = (char) refBase;
                    padm.ALT_BASE = (char) altBase;

                    bbdm.SAMPLE_ALIAS = sampleAlias;
                    bbdm.LIBRARY = library;
                    bbdm.CONTEXT = context;
                    bbdm.REF_BASE = (char) refBase;
                    bbdm.ALT_BASE = (char) altBase;

                    // do the actual counting, as explained in the metric definitions
                    padm.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_NEG + revRefAlignments.R1_NEG + revRefAlignments.R2_POS;
                    padm.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_NEG + revAltAlignments.R1_NEG + revAltAlignments.R2_POS;
                    padm.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + revRefAlignments.R1_POS + revRefAlignments.R2_NEG;
                    padm.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + revAltAlignments.R1_POS + revAltAlignments.R2_NEG;

                    bbdm.FWD_CXT_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + fwdRefAlignments.R2_NEG;
                    bbdm.FWD_CXT_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + fwdAltAlignments.R2_NEG;
                    bbdm.REV_CXT_REF_BASES = revRefAlignments.R1_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_POS + revRefAlignments.R2_NEG;
                    bbdm.REV_CXT_ALT_BASES = revAltAlignments.R1_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_POS + revAltAlignments.R2_NEG;

                    // calculate derived stats (conveniently in a separate method)
                    padm.calculateDerivedStatistics();
                    bbdm.calculateDerivedStatistics();

                    // add to list
                    detailMetrics.preAdaptDetailMetrics.add(padm);
                    detailMetrics.baitBiasDetailMetrics.add(bbdm);
                }
            }
        }
        return detailMetrics;
    }

    /**
     * Core extraction method for summary metrics (which includes extraction of detail metrics).
     */
    private LibraryLevelMetrics extractAllMetrics(final String sampleAlias, final String library, final ContextAccumulator contextSeparatedCounts) {
        // compute detail metrics for each individual reference context
        final DetailMetricsCollector contextLevelMetrics = extractDetailMetrics(sampleAlias, library, contextSeparatedCounts);

        // compute detail metrics for each of the 12 single-base errors
        final DetailMetricsCollector refBaseLevelMetrics;
        final boolean variableContexts = (CONTEXT_SIZE > 0);
        if (variableContexts) {
            final Set<String> contexts = contextSeparatedCounts.keySet();
            final ListMap<String, String> partition = new ListMap<String, String>();
            for (final String context : contexts) {
                final char refBase = getCentralBase(context);
                partition.add(String.valueOf(refBase), context);
            }
            final ContextAccumulator contextCombinedCounts = contextSeparatedCounts.reduceBy(partition);
            refBaseLevelMetrics = extractDetailMetrics(sampleAlias, library, contextCombinedCounts);
        } else {
            refBaseLevelMetrics = contextLevelMetrics;
        }

        // construct the state transition "matrices"
        final Map<Byte, Map<Byte, PreAdaptDetailMetrics>> preAdaptArtifacts = new HashMap<Byte, Map<Byte, PreAdaptDetailMetrics>>();
        final Map<Byte, Map<Byte, BaitBiasDetailMetrics>> baitBiasArtifacts = new HashMap<Byte, Map<Byte, BaitBiasDetailMetrics>>();
        for (final byte b : BASES) {
            preAdaptArtifacts.put(b, new HashMap<Byte, PreAdaptDetailMetrics>());
            baitBiasArtifacts.put(b, new HashMap<Byte, BaitBiasDetailMetrics>());
        }
        for (final PreAdaptDetailMetrics padm : refBaseLevelMetrics.preAdaptDetailMetrics) {
            preAdaptArtifacts.get((byte) padm.REF_BASE).put((byte) padm.ALT_BASE, padm);
        }
        for (final BaitBiasDetailMetrics bbdm : refBaseLevelMetrics.baitBiasDetailMetrics) {
            baitBiasArtifacts.get((byte) bbdm.REF_BASE).put((byte) bbdm.ALT_BASE, bbdm);
        }

        // compute summary metrics. painfully.
        // TODO could this be avoided by doing something clever with the names of the summary metric fields?
        final PreAdaptSummaryMetrics preAdaptSummaryMetrics = new PreAdaptSummaryMetrics();
        preAdaptSummaryMetrics.SAMPLE_ALIAS = sampleAlias;
        preAdaptSummaryMetrics.LIBRARY = library;
        preAdaptSummaryMetrics.A_TO_C_ERROR_RATE = preAdaptArtifacts.get((byte) 'A').get((byte) 'C').ERROR_RATE;
        preAdaptSummaryMetrics.A_TO_G_ERROR_RATE = preAdaptArtifacts.get((byte) 'A').get((byte) 'G').ERROR_RATE;
        preAdaptSummaryMetrics.A_TO_T_ERROR_RATE = preAdaptArtifacts.get((byte) 'A').get((byte) 'T').ERROR_RATE;
        preAdaptSummaryMetrics.C_TO_A_ERROR_RATE = preAdaptArtifacts.get((byte) 'C').get((byte) 'A').ERROR_RATE;
        preAdaptSummaryMetrics.C_TO_G_ERROR_RATE = preAdaptArtifacts.get((byte) 'C').get((byte) 'G').ERROR_RATE;
        preAdaptSummaryMetrics.C_TO_T_ERROR_RATE = preAdaptArtifacts.get((byte) 'C').get((byte) 'T').ERROR_RATE;
        preAdaptSummaryMetrics.G_TO_A_ERROR_RATE = preAdaptArtifacts.get((byte) 'G').get((byte) 'A').ERROR_RATE;
        preAdaptSummaryMetrics.G_TO_C_ERROR_RATE = preAdaptArtifacts.get((byte) 'G').get((byte) 'C').ERROR_RATE;
        preAdaptSummaryMetrics.G_TO_T_ERROR_RATE = preAdaptArtifacts.get((byte) 'G').get((byte) 'T').ERROR_RATE;
        preAdaptSummaryMetrics.T_TO_A_ERROR_RATE = preAdaptArtifacts.get((byte) 'T').get((byte) 'A').ERROR_RATE;
        preAdaptSummaryMetrics.T_TO_C_ERROR_RATE = preAdaptArtifacts.get((byte) 'T').get((byte) 'C').ERROR_RATE;
        preAdaptSummaryMetrics.T_TO_G_ERROR_RATE = preAdaptArtifacts.get((byte) 'T').get((byte) 'G').ERROR_RATE;
        preAdaptSummaryMetrics.A_TO_C_QSCORE = preAdaptArtifacts.get((byte) 'A').get((byte) 'C').QSCORE;
        preAdaptSummaryMetrics.A_TO_G_QSCORE = preAdaptArtifacts.get((byte) 'A').get((byte) 'G').QSCORE;
        preAdaptSummaryMetrics.A_TO_T_QSCORE = preAdaptArtifacts.get((byte) 'A').get((byte) 'T').QSCORE;
        preAdaptSummaryMetrics.C_TO_A_QSCORE = preAdaptArtifacts.get((byte) 'C').get((byte) 'A').QSCORE;
        preAdaptSummaryMetrics.C_TO_G_QSCORE = preAdaptArtifacts.get((byte) 'C').get((byte) 'G').QSCORE;
        preAdaptSummaryMetrics.C_TO_T_QSCORE = preAdaptArtifacts.get((byte) 'C').get((byte) 'T').QSCORE;
        preAdaptSummaryMetrics.G_TO_A_QSCORE = preAdaptArtifacts.get((byte) 'G').get((byte) 'A').QSCORE;
        preAdaptSummaryMetrics.G_TO_C_QSCORE = preAdaptArtifacts.get((byte) 'G').get((byte) 'C').QSCORE;
        preAdaptSummaryMetrics.G_TO_T_QSCORE = preAdaptArtifacts.get((byte) 'G').get((byte) 'T').QSCORE;
        preAdaptSummaryMetrics.T_TO_A_QSCORE = preAdaptArtifacts.get((byte) 'T').get((byte) 'A').QSCORE;
        preAdaptSummaryMetrics.T_TO_C_QSCORE = preAdaptArtifacts.get((byte) 'T').get((byte) 'C').QSCORE;
        preAdaptSummaryMetrics.T_TO_G_QSCORE = preAdaptArtifacts.get((byte) 'T').get((byte) 'G').QSCORE;

        final BaitBiasSummaryMetrics baitBiasSummaryMetrics = new BaitBiasSummaryMetrics();
        baitBiasSummaryMetrics.SAMPLE_ALIAS = sampleAlias;
        baitBiasSummaryMetrics.LIBRARY = library;
        baitBiasSummaryMetrics.A_TO_C_ERROR_RATE = baitBiasArtifacts.get((byte) 'A').get((byte) 'C').ERROR_RATE;
        baitBiasSummaryMetrics.A_TO_G_ERROR_RATE = baitBiasArtifacts.get((byte) 'A').get((byte) 'G').ERROR_RATE;
        baitBiasSummaryMetrics.A_TO_T_ERROR_RATE = baitBiasArtifacts.get((byte) 'A').get((byte) 'T').ERROR_RATE;
        baitBiasSummaryMetrics.C_TO_A_ERROR_RATE = baitBiasArtifacts.get((byte) 'C').get((byte) 'A').ERROR_RATE;
        baitBiasSummaryMetrics.C_TO_G_ERROR_RATE = baitBiasArtifacts.get((byte) 'C').get((byte) 'G').ERROR_RATE;
        baitBiasSummaryMetrics.C_TO_T_ERROR_RATE = baitBiasArtifacts.get((byte) 'C').get((byte) 'T').ERROR_RATE;
        baitBiasSummaryMetrics.G_TO_A_ERROR_RATE = baitBiasArtifacts.get((byte) 'G').get((byte) 'A').ERROR_RATE;
        baitBiasSummaryMetrics.G_TO_C_ERROR_RATE = baitBiasArtifacts.get((byte) 'G').get((byte) 'C').ERROR_RATE;
        baitBiasSummaryMetrics.G_TO_T_ERROR_RATE = baitBiasArtifacts.get((byte) 'G').get((byte) 'T').ERROR_RATE;
        baitBiasSummaryMetrics.T_TO_A_ERROR_RATE = baitBiasArtifacts.get((byte) 'T').get((byte) 'A').ERROR_RATE;
        baitBiasSummaryMetrics.T_TO_C_ERROR_RATE = baitBiasArtifacts.get((byte) 'T').get((byte) 'C').ERROR_RATE;
        baitBiasSummaryMetrics.T_TO_G_ERROR_RATE = baitBiasArtifacts.get((byte) 'T').get((byte) 'G').ERROR_RATE;
        baitBiasSummaryMetrics.A_TO_C_QSCORE = baitBiasArtifacts.get((byte) 'A').get((byte) 'C').QSCORE;
        baitBiasSummaryMetrics.A_TO_G_QSCORE = baitBiasArtifacts.get((byte) 'A').get((byte) 'G').QSCORE;
        baitBiasSummaryMetrics.A_TO_T_QSCORE = baitBiasArtifacts.get((byte) 'A').get((byte) 'T').QSCORE;
        baitBiasSummaryMetrics.C_TO_A_QSCORE = baitBiasArtifacts.get((byte) 'C').get((byte) 'A').QSCORE;
        baitBiasSummaryMetrics.C_TO_G_QSCORE = baitBiasArtifacts.get((byte) 'C').get((byte) 'G').QSCORE;
        baitBiasSummaryMetrics.C_TO_T_QSCORE = baitBiasArtifacts.get((byte) 'C').get((byte) 'T').QSCORE;
        baitBiasSummaryMetrics.G_TO_A_QSCORE = baitBiasArtifacts.get((byte) 'G').get((byte) 'A').QSCORE;
        baitBiasSummaryMetrics.G_TO_C_QSCORE = baitBiasArtifacts.get((byte) 'G').get((byte) 'C').QSCORE;
        baitBiasSummaryMetrics.G_TO_T_QSCORE = baitBiasArtifacts.get((byte) 'G').get((byte) 'T').QSCORE;
        baitBiasSummaryMetrics.T_TO_A_QSCORE = baitBiasArtifacts.get((byte) 'T').get((byte) 'A').QSCORE;
        baitBiasSummaryMetrics.T_TO_C_QSCORE = baitBiasArtifacts.get((byte) 'T').get((byte) 'C').QSCORE;
        baitBiasSummaryMetrics.T_TO_G_QSCORE = baitBiasArtifacts.get((byte) 'T').get((byte) 'G').QSCORE;

        // wrap up
        final List<PreAdaptDetailMetrics> allPreAdaptDetailMetrics = new ArrayList<PreAdaptDetailMetrics>();
        final List<BaitBiasDetailMetrics> allBaitBiasDetailMetrics = new ArrayList<BaitBiasDetailMetrics>();
        if (variableContexts) {
            allPreAdaptDetailMetrics.addAll(contextLevelMetrics.preAdaptDetailMetrics);
            allBaitBiasDetailMetrics.addAll(contextLevelMetrics.baitBiasDetailMetrics);
        }
        allPreAdaptDetailMetrics.addAll(refBaseLevelMetrics.preAdaptDetailMetrics);
        allBaitBiasDetailMetrics.addAll(refBaseLevelMetrics.baitBiasDetailMetrics);
        return new LibraryLevelMetrics(allPreAdaptDetailMetrics, allBaitBiasDetailMetrics, preAdaptSummaryMetrics, baitBiasSummaryMetrics);
    }

    // --------------------------------------------------------------------------------------------
    // Helper classes
    // --------------------------------------------------------------------------------------------

    /**
     * SAM filter for insert size range.
     */
    private static class InsertSizeFilter implements SamRecordFilter {
        final int minInsertSize;
        final int maxInsertSize;

        InsertSizeFilter(final int minInsertSize, final int maxInsertSize) {
            this.minInsertSize = minInsertSize;
            this.maxInsertSize = maxInsertSize;
        }

        @Override
        public boolean filterOut(final SAMRecord rec) {
            // Treat both parameters == 0 as not filtering
            if (minInsertSize == 0 && maxInsertSize == 0) return false;

            if (rec.getReadPairedFlag()) {
                final int ins = Math.abs(rec.getInferredInsertSize());
                return ins < minInsertSize || ins > maxInsertSize;
            }

            // If the read isn't paired and either min or max is specified filter it out
            return minInsertSize != 0 || maxInsertSize != 0;
        }

        @Override
        public boolean filterOut(final SAMRecord r1, final SAMRecord r2) {
            return filterOut(r1) || filterOut(r2);
        }
    }

    /**
     * Main accumulator class
     */
    private class ArtifactCounter {
        private final Set<String> acceptedContexts;
        private final Set<String> acceptedLibraries;
        private final Map<String, ContextAccumulator> libraryMap;

        private ArtifactCounter(final Set<String> libraries, final Set<String> contexts) {
            this.acceptedLibraries = libraries;
            this.acceptedContexts = contexts;
            this.libraryMap = new HashMap<String, ContextAccumulator>();
            for (final String library : libraries) {
                this.libraryMap.put(library, new ContextAccumulator(contexts));
            }
        }

        private void countAlleles(final SamLocusIterator.LocusInfo info, final String refContext) {
            for (final SamLocusIterator.RecordAndOffset rec : info.getRecordAndPositions()) {
                final byte qual;
                final SAMRecord samrec = rec.getRecord();

                if (USE_OQ) {
                    final byte[] oqs = samrec.getOriginalBaseQualities();
                    if (oqs != null) qual = oqs[rec.getOffset()];
                    else qual = rec.getBaseQuality();
                } else {
                    qual = rec.getBaseQuality();
                }

                // Skip if below qual
                if (qual < MINIMUM_QUALITY_SCORE) continue;

                // Skip if context or library is unknown
                final String library = nvl(samrec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
                if (!acceptedLibraries.contains(library)) continue;
                if (!acceptedContexts.contains(refContext)) continue;

                // Count the base
                final byte readBase = rec.getReadBase();
                this.libraryMap.get(library).addRecord(samrec, readBase, refContext);
            }
        }

        private List<LibraryLevelMetrics> finalize(final String sampleAlias) {
            final List<LibraryLevelMetrics> allMetrics = new ArrayList<LibraryLevelMetrics>();
            for (final String library : this.libraryMap.keySet()) {
                final ContextAccumulator accumulator = this.libraryMap.get(library);
                final LibraryLevelMetrics metrics = extractAllMetrics(sampleAlias, library, accumulator);
                allMetrics.add(metrics);
            }
            return allMetrics;
        }
    }

    /**
     * Breaks down alignments by read1/read2 and positive/negative strand.
     */
    private static class AlignmentAccumulator {
        private long R1_POS = 0;
        private long R1_NEG = 0;
        private long R2_POS = 0;
        private long R2_NEG = 0;

        private void addRecord(final SAMRecord rec) {
            final boolean isNegativeStrand = rec.getReadNegativeStrandFlag();
            final boolean isReadTwo = rec.getReadPairedFlag() && rec.getSecondOfPairFlag();
            if (isReadTwo) {
                if (isNegativeStrand) this.R2_NEG++;
                else this.R2_POS++;
            } else {
                if (isNegativeStrand) this.R1_NEG++;
                else this.R1_POS++;
            }
        }

        private static AlignmentAccumulator combine(final Iterable<AlignmentAccumulator> accumulators) {
            final AlignmentAccumulator combined = new AlignmentAccumulator();
            for (final AlignmentAccumulator acc : accumulators) {
                combined.R1_POS += acc.R1_POS;
                combined.R1_NEG += acc.R1_NEG;
                combined.R2_POS += acc.R2_POS;
                combined.R2_NEG += acc.R2_NEG;
            }
            return combined;
        }
    }

    /**
     * One level above AlignmentAccumulator - separate by called base.
     */
    private static class CalledBaseAccumulator extends HashMap<Byte, AlignmentAccumulator> {
        private CalledBaseAccumulator() {
            for (final byte b : BASES) {
                this.put(b, new AlignmentAccumulator());
            }
        }

        private void addRecord(final SAMRecord rec, final byte readBase) {
            this.get(readBase).addRecord(rec);
        }

        /**
         * Combine a set of CalledBaseAccumulators by summing together the AlignmentAccumulators
         * for each key (i.e. base). This assumes that each keySet contains all four bases.
         */
        private static CalledBaseAccumulator combine(final Iterable<CalledBaseAccumulator> accumulators) {
            final ListMap<Byte, AlignmentAccumulator> inverted = new ListMap<Byte, AlignmentAccumulator>();
            for (final CalledBaseAccumulator acc : accumulators) {
                for (final byte b : BASES) {
                    inverted.add(b, acc.get(b));
                }
            }
            final CalledBaseAccumulator combined = new CalledBaseAccumulator();
            for (final byte b : BASES) {
                combined.put(b, AlignmentAccumulator.combine(inverted.get(b)));
            }
            return combined;
        }
    }

    /**
     * One level above CalledBaseAccumulator - separate by reference context.
     */
    private static class ContextAccumulator extends HashMap<String, CalledBaseAccumulator> {
        private ContextAccumulator(final Iterable<String> allowedContexts) {
            for (final String context : allowedContexts) {
                this.put(context, new CalledBaseAccumulator());
            }
        }

        private void addRecord(final SAMRecord rec, final byte readBase, final String context) {
            this.get(context).addRecord(rec, readBase);
        }

        /**
         * Given a partitioning of the contexts, create a new ContextAccumulator with the contexts' contents
         * summed together according to the partition (whew!)
         *
         * Note: if the input refers to a context that isn't found in the current ContextAccumulator, we simply ignore it.
         *
         */
        private ContextAccumulator reduceBy(final ListMap<String, String> partition) {
            final Set<String> keys = partition.keySet();
            final ContextAccumulator combined = new ContextAccumulator(keys);
            for (final String key : keys) {
                final List<CalledBaseAccumulator> toBeCombined = new ArrayList<CalledBaseAccumulator>();
                for (final String context : partition.get(key)) {
                    if (this.containsKey(context)) {
                        toBeCombined.add(this.get(context));
                    }
                }
                combined.put(key, CalledBaseAccumulator.combine(toBeCombined));
            }
            return combined;
        }
    }

    /** A container for metrics at the per-library level. */
    private static class LibraryLevelMetrics {
        private final List<PreAdaptDetailMetrics> preAdaptDetailMetrics;
        private final List<BaitBiasDetailMetrics> baitBiasDetailMetrics;
        private final PreAdaptSummaryMetrics preAdaptSummaryMetrics;
        private final BaitBiasSummaryMetrics baitBiasSummaryMetrics;

        private LibraryLevelMetrics(final List<PreAdaptDetailMetrics> preAdaptDetailMetrics,
                                    final List<BaitBiasDetailMetrics> baitBiasDetailMetrics,
                                    final PreAdaptSummaryMetrics preAdaptSummaryMetrics,
                                    final BaitBiasSummaryMetrics baitBiasSummaryMetrics) {
            this.preAdaptDetailMetrics = preAdaptDetailMetrics;
            this.baitBiasDetailMetrics = baitBiasDetailMetrics;
            this.preAdaptSummaryMetrics = preAdaptSummaryMetrics;
            this.baitBiasSummaryMetrics = baitBiasSummaryMetrics;
        }
    }

    private static class DetailMetricsCollector {
        private final List<PreAdaptDetailMetrics> preAdaptDetailMetrics = new ArrayList<PreAdaptDetailMetrics>();
        private final List<BaitBiasDetailMetrics> baitBiasDetailMetrics = new ArrayList<BaitBiasDetailMetrics>();
    }

    // --------------------------------------------------------------------------------------------
    // Helper methods
    // --------------------------------------------------------------------------------------------

    /** Mimic of Oracle's nvl() - returns the first value if not null, otherwise the second value. */
    private <T> T nvl(final T value1, final T value2) {
        if (value1 != null) return value1;
        else return value2;
    }

    /**
     * Little method to expand a set of sequences to include all reverse complements.
     * Necessary if a user manually specifies a set of contexts, due to symmetries that the analysis depends on.
     */
    private Set<String> includeReverseComplements(final Set<String> sequences) {
        final Set<String> all = new HashSet<String>();
        for (final String seq : sequences) {
            all.add(seq);
            all.add(SequenceUtil.reverseComplement(seq));
        }
        return all;
    }

    private Set<String> makeContextStrings(final int length) {
        final Set<String> contexts = new HashSet<String>();

        for (final byte[] kmer : generateAllKmers(length)) {
            contexts.add(StringUtil.bytesToString(kmer));
        }

        log.info("Generated " + contexts.size() + " context strings.");
        return contexts;
    }

    /** Generates all possible unambiguous kmers of length and returns them as byte[]s. */
    private List<byte[]> generateAllKmers(final int length) {
        final List<byte[]> sofar = new LinkedList<byte[]>();

        if (sofar.size() == 0) {
            sofar.add(new byte[length]);
        }

        while (true) {
            final byte[] bs = sofar.remove(0);
            final int indexOfNextBase = findIndexOfNextBase(bs);

            if (indexOfNextBase == -1) {
                sofar.add(bs);
                break;
            } else {
                for (final byte b : BASES) {
                    final byte[] next = Arrays.copyOf(bs, bs.length);
                    next[indexOfNextBase] = b;
                    sofar.add(next);
                }
            }
        }

        return sofar;
    }

    /** Finds the first zero character in the array, or returns -1 if all are non-zero. */
    private int findIndexOfNextBase(final byte[] bs) {
        for (int i = 0; i < bs.length; ++i) {
            if (bs[i] == 0) return i;
        }

        return -1;
    }

    /** Get the central base of a context string. Throw an error if the string length is even. */
    private char getCentralBase(final String context) {
        if (context.length() % 2 == 0) throw new PicardException("Context " + context + " has an even length, and thus no central base.");
        else return context.charAt(context.length() / 2);
    }

    // --------------------------------------------------------------------------------------------
    // Metrics classes
    // --------------------------------------------------------------------------------------------

    /**
     * Library-level summary of estimated DNA damage occurring before adapters are added.
     *
     * These errors occur on the original template strand, and thus correlate with
     * read number / orientation in a specific way.
     */
    public static class PreAdaptSummaryMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;

        /** A>C error rate. */
        public double A_TO_C_ERROR_RATE;
        /** A>C Phred-scaled quality score. */
        public double A_TO_C_QSCORE;
        /** A>G error rate. */
        public double A_TO_G_ERROR_RATE;
        /** A>G Phred-scaled quality score. */
        public double A_TO_G_QSCORE;
        /** A>T error rate. */
        public double A_TO_T_ERROR_RATE;
        /** A>T Phred-scaled quality score. */
        public double A_TO_T_QSCORE;

        /** C>A error rate. */
        public double C_TO_A_ERROR_RATE;
        /** C>A Phred-scaled quality score. */
        public double C_TO_A_QSCORE;
        /** C>G error rate. */
        public double C_TO_G_ERROR_RATE;
        /** C>G Phred-scaled quality score. */
        public double C_TO_G_QSCORE;
        /** C>T error rate. */
        public double C_TO_T_ERROR_RATE;
        /** C>T Phred-scaled quality score. */
        public double C_TO_T_QSCORE;

        /** G>A error rate. */
        public double G_TO_A_ERROR_RATE;
        /** G>A Phred-scaled quality score. */
        public double G_TO_A_QSCORE;
        /** G>C error rate. */
        public double G_TO_C_ERROR_RATE;
        /** G>C Phred-scaled quality score. */
        public double G_TO_C_QSCORE;
        /** G>T error rate. */
        public double G_TO_T_ERROR_RATE;
        /** G>T Phred-scaled quality score. */
        public double G_TO_T_QSCORE;

        /** T>A error rate. */
        public double T_TO_A_ERROR_RATE;
        /** T>A Phred-scaled quality score. */
        public double T_TO_A_QSCORE;
        /** T>C error rate. */
        public double T_TO_C_ERROR_RATE;
        /** T>C Phred-scaled quality score. */
        public double T_TO_C_QSCORE;
        /** T>G error rate. */
        public double T_TO_G_ERROR_RATE;
        /** T>G Phred-scaled quality score. */
        public double T_TO_G_QSCORE;
    }

    /**
     * Library-level summary of estimated DNA damage occurring after baits are added.
     *
     * These occur on the strand being targeted by the bait, and correlate with some
     * kind of substitution being observed preferentially on the strand relative to
     * its complement (e.g. G->T vs C->A).
     */
    public static class BaitBiasSummaryMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;

        /** A>C error rate. */
        public double A_TO_C_ERROR_RATE;
        /** A>C Phred-scaled quality score. */
        public double A_TO_C_QSCORE;
        /** A>G error rate. */
        public double A_TO_G_ERROR_RATE;
        /** A>G Phred-scaled quality score. */
        public double A_TO_G_QSCORE;
        /** A>T error rate. */
        public double A_TO_T_ERROR_RATE;
        /** A>T Phred-scaled quality score. */
        public double A_TO_T_QSCORE;

        /** C>A error rate. */
        public double C_TO_A_ERROR_RATE;
        /** C>A Phred-scaled quality score. */
        public double C_TO_A_QSCORE;
        /** C>G error rate. */
        public double C_TO_G_ERROR_RATE;
        /** C>G Phred-scaled quality score. */
        public double C_TO_G_QSCORE;
        /** C>T error rate. */
        public double C_TO_T_ERROR_RATE;
        /** C>T Phred-scaled quality score. */
        public double C_TO_T_QSCORE;

        /** G>A error rate. */
        public double G_TO_A_ERROR_RATE;
        /** G>A Phred-scaled quality score. */
        public double G_TO_A_QSCORE;
        /** G>C error rate. */
        public double G_TO_C_ERROR_RATE;
        /** G>C Phred-scaled quality score. */
        public double G_TO_C_QSCORE;
        /** G>T error rate. */
        public double G_TO_T_ERROR_RATE;
        /** G>T Phred-scaled quality score. */
        public double G_TO_T_QSCORE;

        /** T>A error rate. */
        public double T_TO_A_ERROR_RATE;
        /** T>A Phred-scaled quality score. */
        public double T_TO_A_QSCORE;
        /** T>C error rate. */
        public double T_TO_C_ERROR_RATE;
        /** T>C Phred-scaled quality score. */
        public double T_TO_C_QSCORE;
        /** T>G error rate. */
        public double T_TO_G_ERROR_RATE;
        /** T>G Phred-scaled quality score. */
        public double T_TO_G_QSCORE;
    }

    /**
     * Pre-adapter DNA damage broken down by context (the reference bases surrounding the base of interest).
     */
    public static class PreAdaptDetailMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The sequence context to which the analysis is constrained. */
        public String CONTEXT;

        /** The original base on the reference strand. */
        public char REF_BASE;
        /** The alternative base that is "inferred" as a result of DNA damage. */
        public char ALT_BASE;

        /** The number of REF_BASE:REF_BASE alignments having a read number and orientation that supports the presence of this artifact. */
        public long PRO_REF_BASES;
        /** The number of REF_BASE:ALT_BASE alignments having a read number and orientation that supports the presence of this artifact. */
        public long PRO_ALT_BASES;
        /** The number of REF_BASE:REF_BASE alignments having a read number and orientation that refutes the presence of this artifact. */
        public long CON_REF_BASES;
        /** The number of REF_BASE:ALT_BASE alignments having a read number and orientation that refutes the presence of this artifact. */
        public long CON_ALT_BASES;

        /**
         * The estimated error rate due to this artifact.
         * Calculated as max(1e-10, (PRO_ALT_BASES - CON_ALT_BASES) / (PRO_ALT_BASES + PRO_REF_BASES + CON_ALT_BASES + CON_REF_BASES)).
         */
        public double ERROR_RATE;
        /** The Phred-scaled quality score of this artifact, calculated as -10 * log10(ERROR_RATE). */
        public double QSCORE;

        /**
         * Calculate the error rate given the raw counts. Negative rates are set to MIN_ERROR.
         *
         * TODO explain this calculation
         */
        private void calculateDerivedStatistics() {
            this.ERROR_RATE = MIN_ERROR;
            final long totalBases = this.PRO_REF_BASES + this.PRO_ALT_BASES + this.CON_REF_BASES + this.CON_ALT_BASES;
            if (totalBases > 0) {
                final double rawErrorRate = (this.PRO_ALT_BASES - this.CON_ALT_BASES) / (double) totalBases;
                this.ERROR_RATE = Math.max(MIN_ERROR, rawErrorRate);
            }
            this.QSCORE = -10 * log10(this.ERROR_RATE);
        }
    }

    /**
     * Bait bias DNA damage broken down by context (the reference bases surrounding the base of interest).
     */
    public static class BaitBiasDetailMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The sequence context to which the analysis is constrained. */
        public String CONTEXT;

        /** The reference base being reported on. */
        public char REF_BASE;
        /** The alternative base being reported on. */
        public char ALT_BASE;

        /** The number of REF_BASE:REF_BASE alignments at sites with the given reference context. */
        public long FWD_CXT_REF_BASES;
        /** The number of REF_BASE:ALT_BASE alignments at sites with the given reference context. */
        public long FWD_CXT_ALT_BASES;
        /** The number of ~REF_BASE:~REF_BASE alignments at sites complementary to the given reference context. */
        public long REV_CXT_REF_BASES;
        /** The number of ~REF_BASE:~ALT_BASE alignments at sites complementary to the given reference context. */
        public long REV_CXT_ALT_BASES;

        /** The error rate of REF_BASE:ALT_BASE, calculated as max(1e-10, FWD_CXT_ALT_BASES / (FWD_CXT_ALT_BASES + FWD_CXT_REF_BASES)). */
        public double FWD_ERROR_RATE;
        /** The error rate of ~REF_BASE:~ALT_BASE, calculated as max(1e-10, REV_CXT_ALT_BASES / (REV_CXT_ALT_BASES + REV_CXT_REF_BASES)). */
        public double REV_ERROR_RATE;

        /** The combined error rate (a.k.a. reference base bias), calculated as max(1e-10, FWD_ERROR_RATE - REV_ERROR_RATE). */
        public double ERROR_RATE;
        /** The Phred-scaled quality score of the reference base bias, calculated as -10 * log10(ERROR_RATE). */
        public double QSCORE;

        /**
         * Calculate the error rate given the raw counts. Negative rates are set to MIN_ERROR.
         *
         * TODO explain this calculation
         */
        private void calculateDerivedStatistics() {

            // TODO are the denominators in these error rates correct?

            this.FWD_ERROR_RATE = MIN_ERROR;
            final long fwdBases = this.FWD_CXT_REF_BASES + this.FWD_CXT_ALT_BASES;
            if (fwdBases > 0) {
                final double fwdErr = this.FWD_CXT_ALT_BASES / (double) fwdBases;
                this.FWD_ERROR_RATE = Math.max(MIN_ERROR, fwdErr);
            }

            this.REV_ERROR_RATE = MIN_ERROR;
            final long revBases = this.REV_CXT_REF_BASES + this.REV_CXT_ALT_BASES;
            if (revBases > 0) {
                final double revErr = this.REV_CXT_ALT_BASES / (double) revBases;
                this.REV_ERROR_RATE = Math.max(MIN_ERROR, revErr);
            }

            this.ERROR_RATE = Math.max(MIN_ERROR, this.FWD_ERROR_RATE - this.REV_ERROR_RATE);
            this.QSCORE = -10 * log10(this.ERROR_RATE);
        }
    }
}
