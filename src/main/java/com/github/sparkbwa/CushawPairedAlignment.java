/**
 * Copyright 2016 José Manuel Abuín Mosquera <josemanuel.abuin@usc.es>
 *
 * <p>This file is part of SparkBWA.
 *
 * <p>SparkBWA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * <p>SparkBWA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * <p>You should have received a copy of the GNU General Public License along with SparkBWA. If not,
 * see <http://www.gnu.org/licenses/>.
 */

package com.github.sparkbwa;

import org.apache.spark.SparkContext;
import org.apache.spark.api.java.function.Function2;
import scala.Tuple2;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * Class to perform the alignment over a split from the RDD of paired reads using Cushaw as aligner
 *
 * @author José M. Abuín
 */
public class CushawPairedAlignment extends CushawAlignmentBase implements Function2<Integer, Iterator<Tuple2<String, String>>, Iterator<String>> {

    /**
     * Constructor
     * @param context The Spark context
     * @param indexPath Path to the index to use for alignment
     */
    public CushawPairedAlignment(SparkContext context, String indexPath) {
        super(context, indexPath);
        this.cushaw = new CushawJni();
    }

    /**
     * Code to run in each one of the mappers. This is, the alignment with the corresponding entry
     * data The entry data has to be written into the local filesystem
     * @param arg0 The RDD Id
     * @param arg1 An iterator containing the values in this RDD
     * @return An iterator containing the sam file name generated
     * @throws Exception
     */
    public Iterator<String> call(Integer arg0, Iterator<Tuple2<String, String>> arg1) throws Exception {

        LOG.info("["+this.getClass().getName()+"] :: Tmp dir: " + this.tmpDir);

        String inputSequences1 = "";
        String inputSequences2 = "";

        ArrayList<String> sequences1 = new ArrayList<String>();
        ArrayList<String> sequences2 = new ArrayList<String>();

        ArrayList<String> results = new ArrayList<String>();

        String tmp1;
        String tmp2;
        //Iterator<Tuple2<String, String>> sequences = arg1;

        while (arg1.hasNext()) {
            tmp1 = arg1.next()._1;
            tmp2 = arg1.next()._2;

            inputSequences1 = inputSequences1 + tmp1;
            inputSequences2 = inputSequences2 + tmp2;
            sequences1.add(tmp1);
            sequences2.add(tmp2);
        }

        this.cushaw.executeEstimateJava(inputSequences1, inputSequences2);

        for(int i = 0; i< sequences1.size(); i++) {

            results.add(this.cushaw.alignJava(sequences1.get(i), sequences2.get(2)));

        }

        return results.iterator();

    }


}
