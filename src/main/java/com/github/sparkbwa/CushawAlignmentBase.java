package com.github.sparkbwa;
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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.SparkContext;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * Abstract class that represents an alignment
 * @author Jose M. Abuin
 */
public abstract class CushawAlignmentBase implements Serializable {

    protected static final Log LOG = LogFactory.getLog(CushawAlignmentBase.class);

    protected String appName	= "";
    protected String appId		= "";
    protected String tmpDir		= "";
    protected String indexPath  = "";
    CushawJni cushaw;
    /**
     * Constructor for this class
     *
     * @brief This constructor creates a BwaAlignment object to process in each one of the mappers
     * @param context The SparkContext to use
     * @param bwaInterpreter The Bwa object used to perform the alignment
     */
    public CushawAlignmentBase(SparkContext context, String indexPath) {

        this.appId			= context.applicationId();
        this.appName		= context.appName();
        this.tmpDir			= context.getLocalProperty("spark.local.dir");
        this.indexPath      = indexPath;


        this.LOG.info("["+this.getClass().getName()+"] :: " + this.appId + " - " + this.appName);
    }

    /**
     *
     * @param readBatchID Identification for the sam file
     * @return A String for the sam file name
     */
    public String getOutputSamFilename(Integer readBatchID) {
        return this.appName + "-" + this.appId + "-" + readBatchID + ".sam";
    }

}
