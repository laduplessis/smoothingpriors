package beast.core.util;

import beast.core.*;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.IllegalFormatCodePointException;
import java.util.List;

/**
 *
 * Always return times from present to the past (this should be relaxed in future, but then the origin should be passed to it)
 *
 * MRCA - present
 * Oldest sample - present
 *
 * Specific dates
 *
 * Branching events (TODO)
 * Sampling events (TODO)
 * Branching and sampling events (TODO)
 * If nr intervals are provided adds an equal number of events into each interval (TODO - almost done, not tested)
 *
 * includeFirst means that the first interval (eg. from root to MRCA is merged with the next).
 *
 * stopAtInput
 * breakAtInput
 * includeLast
 *
 * values input needs to be specified (because of inheritance), but does nothing. (TODO can we make this more elegant?)
 *
 *
 *
 * Created by dlouis on 04/04/16.
 */
public class TreeSlicer extends RealParameter {

    final static int MRCA          = 0,
                     LASTSAMPLE    = 1,
                     EQUIDISTANT   = 2,
                     DATES         = 3,
                     HEIGHTS       = 4;


    final double eps = 1e-6;

    public Input<Tree> treeInput =
        new Input<>("tree", "Tree over which to calculate the slice", Input.Validate.REQUIRED);

    public Input<String> stopInput =
        new Input<>("stop", "Breakpoint to stop the slicing intervals (tmrca/oldestsample/present)", "tmrca");

    public Input<String> typeInput =
            new Input<>("type", "Type of slice (equidistant/dates/heights)", "equidistant");

    // Need to add includeStart
    public Input<Boolean> includeLastInput =
            new Input<>("includeLast", "Include the last breakpoint (stopping criterion) into the vector", true);


    protected Tree tree;
    protected int stop, type, start;        // Where to stop and start the slice
    protected boolean includeLast;
    protected Double [] dates;

    protected double tmrca,                 // Height of TMRCA of the tree
                     oldest,                // Height of the oldest sample
                     newest,                // Height of the most recent sample
                     anchordate;            // Date for translating height to calendar date (most recent sample at time = 0)

    private boolean timesKnown;


    @Override
    public void initAndValidate() {
        int dimension;
        String stopStr, typeStr;


        /* Read tree and anchor date */
        tree = treeInput.get();

        double minheight = tree.getRoot().getHeight();
        for (Node N : tree.getNodesAsArray()) {
            if (N.getHeight() <= minheight) {
                anchordate = N.getDate();
                minheight  = N.getHeight();
            }
        }


        /* Read inputs */
        if (stopInput.get() != null)
            stopStr = stopInput.get().toLowerCase().trim();
        else
            stopStr = "tmrca";

        if (typeInput.get() != null)
            typeStr = typeInput.get().toLowerCase().trim();
        else
            typeStr = "equidistant";


        /* Initialize values and set type */
        Double[] valuesRaw = valuesInput.get().toArray(new Double[0]);

        if (typeStr.equals("dates")) {
            type = DATES;

            dates = new Double [valuesRaw.length];
            for (int i = 0; i < valuesRaw.length; i++)
                dates[i] = valuesRaw[i];
            HeapSort.sort(dates);

            if (dates[dates.length-1] == anchordate) {
                dimension = dates.length;
            } else
                dimension = dates.length+1;

            dimensionInput.setValue(dimension, this);
        } else
        if (typeStr.equals("equidistant")) {
            type = EQUIDISTANT;
            dimension = dimensionInput.get();
        } else
            throw new IllegalArgumentException("Unknown treeslice type");

        values = new Double[dimension];
        storedValues = new Double[dimension];



        /* Set stop criterion */
        if (stopStr.equals("tmrca")) {
            stop = MRCA;
        } else
        if (stopStr.equals("oldestsample")) {
            stop = LASTSAMPLE;
        } else
            throw new IllegalArgumentException("Unknown stop criterion!");


        if (includeLastInput.get() != null)
            includeLast = includeLastInput.get();


        calculateTimes(tree);

        // Initialization accounting (not really used)
        m_fLower = Double.NEGATIVE_INFINITY;
        m_fUpper = Double.POSITIVE_INFINITY;
        m_bIsDirty = new boolean[dimensionInput.get()];
        minorDimension = minorDimensionInput.get();
        if (minorDimension > 0 && dimensionInput.get() % minorDimension > 0) {
            throw new IllegalArgumentException("Dimension must be divisible by stride");
        }
        this.storedValues = values.clone();

        timesKnown = false;
    }



    protected double dateToHeight(double date) {
        return (anchordate - date);
    }

    protected double heightToDate(double height) {
        return (anchordate - height);
    }


    /**
     *
     * Update a few state variables
     *
     *  - tmrca   (height of tmrca if present = 0 ... tmrca_date = present - tmrca
     *  - oldest  (height of the oldest sample)
     *  - newest  (height of the youngest sample - should be 0)
     *
     * O(n) for n nodes
     *
     * @param tree
     */
    protected void updateAnchorTimes(Tree tree) {
        double height;

        tmrca = tree.getRoot().getHeight();

        oldest = 0;
        newest = tmrca;
        for (Node N : tree.getNodesAsArray()) {
            if (N.isLeaf()) {
                height = N.getHeight();

                if (height > oldest) {
                    oldest = height;
                }

                if (height  < newest) {
                    newest  = height;
                }
            }
        }
    }


    // Anchor time does not need to be recalculated if tip-dating is used since the tree is anchored based on
    // the sampling times in the TraitSet (at startup). If tip-dating is used nodes may have negative tip-dates, but the
    // translation between height and date does not change.
    //
    // Sampling times need to be recalculated always, since if tip-dating is used the time of the oldest sample can change
    protected void calculateTimes(Tree tree) {
        double endtime, step;

        updateAnchorTimes(tree);
        //System.out.println(tmrca+"\t"+oldest+"\t"+heightToDate(oldest)+"\t"+newest+"\t"+anchordate);

        if (type == DATES) {

            for (int i = 1; i < getDimension(); i++)
               values[values.length-i] = dateToHeight(dates[i-1]);
            values[0] = 0.0;


        } else {

            if (stop == MRCA) {
                endtime = tmrca;
            } else if (stop == LASTSAMPLE) {
                endtime = oldest + eps;
            } else
                endtime = 1.0;


            if (includeLast)
                step = endtime / (getDimension() - 1);
            else
                step = endtime / (getDimension());


            for (int i = 0; i < getDimension(); i++) {
                values[i] = i * step;
            }

        }



            //System.out.println("MRCA: "+mrca+"\tStepsize: "+step);

//        TraitSet dates = tree.getDateTrait();
//        String[] treetaxa = tree.getTaxaNames();
//
//        for (Node n : tree.getNodesAsArray())
//            if (n.isLeaf())
//                System.out.println(n.getNr()+"\t"+n.getHeight()+'\t'+
//                        treetaxa[n.getNr()]+'\t'+dates.getValue(treetaxa[n.getNr()]));
/*

        for (double value : values)
            System.out.print(heightToDate(value)+"\t");
        System.out.println();

        for (double value : values)
            System.out.print(value+"\t");
        System.out.println();
        System.out.println(this.getDimension());
*/

        timesKnown = true;
    }


    protected boolean requiresRecalculation() {
        timesKnown = false;
        return true;
        //return tree.isDirtyCalculation();
    }


    // Override methods to make sure times get recalculated!

    @Override
    public Double getValue() {
        if (!timesKnown) {
            calculateTimes(tree);
        }
        return values[0];
    }

    @Override
    public Double getValue(final int index) {
        if (!timesKnown) {
            calculateTimes(tree);
        }
        return values[index];
    }

    @Override
    public double getArrayValue() {
        if (!timesKnown) {
            calculateTimes(tree);
        }
        return values[0];
    }

    @Override
    public double getArrayValue(final int index) {
        if (!timesKnown) {
            calculateTimes(tree);
        }
        return values[index];
    }

    @Override
    public Double [] getValues() {
        if (!timesKnown) {
            calculateTimes(tree);
        }
        return Arrays.copyOf(values, values.length);
    }


}
