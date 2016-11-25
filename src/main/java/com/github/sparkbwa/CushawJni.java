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

import cz.adamh.utils.NativeUtils;

import java.io.IOException;

/**
 * Class that calls Cushaw functions by means of JNI
 *
 * @author José M. Abuín
 */
public class CushawJni {

    static {
        try {
            NativeUtils.loadLibraryFromJar("/libcushaw.so");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    //Declaration of native method
    private static native int CushawInit(int argc, String[] argv, int[] lenStrings);
    private static native int parseSequence1(String seq1);
    private static native int parseSequence2(String seq1);
}
