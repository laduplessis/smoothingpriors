package beast.core.util;

import beast.core.*;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;

import java.util.Arrays;

/**
 *
 * Always return times from present to the past
 *
 * MRCA - present
 * Oldest sample - present
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
                     BRANCHES      = 2,
                     SAMPLES       = 3,
                     BRANCHSAMPLES = 4;

    final double eps = 1e-6;


    public Input<Tree> treeInput =
        new Input<>("tree", "Tree over which to calculate the slice", Input.Validate.REQUIRED);

    public Input<String> stopAtInput =
        new Input<>("stopAt", "Where to stop the intervals (mrca/lastsample)", "mrca");

    // Not implemented
    public Input<String> breakAtInput =
            new Input<>("breakAt", "Where to break the intervals (branches/samples/branchsamples)");

    public Input<Boolean> includeLastInput =
            new Input<>("includeLast", "Include the last breakpoint (stopping criterion) into the vector", true);


    protected Tree tree;
    protected int stopCriterion,
                  breakCriterion = -1;
    protected boolean includeLast;
    protected double maxSample;

    private boolean timesKnown;


    @Override
    public void initAndValidate() {
        String stopAt, breakAt;

        super.initAndValidate();

        tree = treeInput.get();

        // Set stop criterion
        if (stopAtInput.get() != null) {
            stopAt = stopAtInput.get().toLowerCase().trim();

            if (stopAt.equals("mrca")) {
                stopCriterion = MRCA;
            } else
            if (stopAt.equals("lastsample")) {
                stopCriterion = LASTSAMPLE;

                // Since the sampling times are fixed we don't need to recalculate them if the tree is dirty
                TraitSet dates = tree.getDateTrait();
                String[] treetaxa = tree.getTaxaNames();

                maxSample = 0;
                for (String taxon : treetaxa) {
                    if (dates.getValue(taxon) > maxSample)
                        maxSample = dates.getValue(taxon);
                }

            } else
                throw new IllegalArgumentException("Unknown stop criterion!");
        }

        // Set break points
        if (breakAtInput.get() != null) {
            breakAt = breakAtInput.get().toLowerCase().trim();

            if (breakAt.equals("branches") || breakAt.equals("branch")) {
                breakCriterion = BRANCHES;
            } else
            if (breakAt.equals("samples") || breakAt.equals("sample")) {
                breakCriterion = SAMPLES;
            } else
            if (breakAt.equals("branchsamples") || breakAt.equals("branchsample")) {
                breakCriterion = BRANCHSAMPLES;
            } else
                throw new IllegalArgumentException("Unknown break criterion!");
        }

        if (includeLastInput.get() != null)
            includeLast = includeLastInput.get();

        calculateTimes(tree);
        timesKnown = false;
    }



    protected void calculateTimes(Tree tree) {
        double endtime, step;


        if (stopCriterion == MRCA) {
            endtime = tree.getRoot().getHeight();
        } else
        if (stopCriterion == LASTSAMPLE) {
            endtime = maxSample+eps;
        } else
            endtime = 1.0;


        if (breakCriterion == BRANCHES) {

            double [] times = new double[tree.getNodeCount()];
            double time;
            int branchevents, groupsize, i, n;

            if (includeLast)
                n = getDimension()-1;
            else
                n = getDimension();

            // Get branching times
            i = 0;
            for (Node N : tree.getNodesAsArray()) {
                time = N.getHeight();
                times[i] = time <= endtime ? time : 0;
                //times[i] = n.getHeight();
                i++;
            }
            HeapSort.sort(times);



            // Count how many branch events between endtime and present
            i = 0;
            while (times[i] <= 0) i++;
            branchevents = times.length-i;

            groupsize = (int) Math.floor((double) branchevents/n);

            values[0] = 0.0;
            for (int j = 1; j < n; j++) {
                values[j] = (times[i]+times[i+1])/2;
                i += groupsize;
            }

            if (includeLast) values[n] = endtime;

            //System.out.print(String.format("1: %.2f", endtime));
            //for (int j = 0; j < getDimension(); j++) {
            //    System.out.print("\t"+String.format("%.2f",values[j]));
            //}
            //System.out.println();
        } else {

            if (includeLast)
                step = endtime/(getDimension()-1);
            else
                step = endtime/(getDimension());

            for (int i = 0; i < getDimension(); i++)
                values[i] = i*step;
        }





            //System.out.println("MRCA: "+mrca+"\tStepsize: "+step);

//        TraitSet dates = tree.getDateTrait();
//        String[] treetaxa = tree.getTaxaNames();
//
//        for (Node n : tree.getNodesAsArray())
//            if (n.isLeaf())
//                System.out.println(n.getNr()+"\t"+n.getHeight()+'\t'+
//                        treetaxa[n.getNr()]+'\t'+dates.getValue(treetaxa[n.getNr()]));


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
