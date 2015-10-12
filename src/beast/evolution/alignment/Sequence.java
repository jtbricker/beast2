/*
* File Sequence.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package beast.evolution.alignment;

import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.util.HashMap;

import beast.core.Description;
import beast.core.Input;
import beast.core.BEASTObject;
import beast.core.Citation;
import beast.evolution.datatype.DataType;

@Citation("Bricker, Justin.  Master's Thesis.  Florida State University, Department of Scientific Computing")
@Description("Single sequence in an alignment WITH QUALITY DATA from FASTQ.")
public class Sequence extends BEASTObject {
    public Input<Integer> totalCountInput = new Input<Integer>("totalcount", "number of states or the number of lineages for this species in SNAPP analysis");
    public Input<String> taxonInput = new Input<String>("taxon", "name of this species", Input.Validate.REQUIRED);
    public Input<String> dataInput = new Input<String>("sequenceData",
            "sequence data, either encoded as a string or as comma separated list of integers, or comma separated likelihoods/probabilities for each site if uncertain=true." +
            "In either case, whitespace is ignored.", Input.Validate.REQUIRED);
    public Input<String> qualityInput = new Input<String>("value",
            "sequence data, either encoded as a string or as comma separated list of integers, or comma separated likelihoods/probabilities for each site if uncertain=true." +
            "In either case, whitespace is ignored.");
    

    
    protected double[][] likelihoods = null;    
    public double[][] getLikelihoods() {
    	return likelihoods;
    }
    
    public Sequence() {
    }

    /**
     * Constructor for testing.
     *
     * @param taxon
     * @param sequence
     * @throws Exception
     */
    public Sequence(String taxon, String sequence) throws Exception {
        taxonInput.setValue(taxon, this);
        dataInput.setValue(sequence, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() throws Exception {
    	initProbabilities();    		
    } // initAndValidate
    
    public void initProbabilities() throws Exception {
    	   	
    	String sdata = dataInput.get();
    	String qdata = qualityInput.get();
    	
        // remove spaces
        sdata = sdata.replaceAll("\\s", "");
        String sStr = sdata.trim();
        char[] sequence = sStr.toCharArray();
        
      //ACGT
  		HashMap<Character, Integer> bases = new HashMap<Character, Integer>();
  	    bases.put('A', 0);
  	    bases.put('C', 1);
  	    bases.put('G', 2);
  	    bases.put('T', 3);
        
        if(qdata != null){
        	qdata = qdata.replaceAll("\\s", "");
        	String qStr = qdata.trim();
        	char[] quality = qStr.toCharArray();
        	
        	double error;
    		for(int i=0; i<sequence.length; i++){
    			error = qualToLikelihood(quality[i]);
    			for(int j=0; j<4; j++){
    				if (likelihoods == null) likelihoods = new double[sequence.length][4];
    				likelihoods[i][j] = error/3;
    			}
    			int seqCharIndex = bases.get(sequence[i]);
    			likelihoods[i][seqCharIndex] = 1-error;
    		}
        }
        else{
        	for(int i=0; i<sequence.length; i++){
    			for(int j=0; j<4; j++){
    				if (likelihoods == null) likelihoods = new double[sequence.length][4];
    				likelihoods[i][j] = 0;
    			}
    			int seqCharIndex = bases.get(sequence[i]);
    			likelihoods[i][seqCharIndex] = 1;
        	}
        }
		
		
		
		
    }

    private double qualToLikelihood(char c) {
    	double x = Math.pow(10,-(double)((int) c - 33)/10);
    	//System.out.println(x);
		return x;
	}

	public List<Integer> getSequence(DataType dataType) throws Exception {
        
    	List<Integer> sequence;
    	if (true) { //true=uncertain
            sequence = new ArrayList<Integer>();
            for (int i=0; i<likelihoods.length; i++) {
            	double m = likelihoods[i][0];
            	int index = 0;
            	for (int j=0; j<likelihoods[i].length; j++) {
            		if (likelihoods[i][j] > m ) {
            			m = likelihoods[i][j];
            			index = j;
            		}        		
            	}
            	sequence.add(index);
            }
    	}

        if (totalCountInput.get() == null) {
            // derive default from char-map
            totalCountInput.setValue(dataType.getStateCount(), this);
        }
        return sequence;
    }

    /**
     * @return the taxon of this sequence as a string.
     */
    public final String getTaxon() {
        return taxonInput.get();
    }

    /**
     * @return the data of this sequence as a string.
     */
    public final String getData() {
        return dataInput.get();
    }
    
    
    /**
     * @return the quality data of this sequence as a string.
     */
    public final String getQuality() {
        return qualityInput.get();
    }


    int mapCharToData(String dataMap, char c) {
        int i = dataMap.indexOf(c);
        if (i >= 0) {
            return i;
        }
        return dataMap.length();
    } // mapCharToData

    /**
     * @param id of target sequence
     * @param sequences a collection of sequences
     * @return the sequence in the collection with the given ID, or null if its not in the collection.
     */
    public static Sequence getSequenceByTaxon(String id, Collection<Sequence> sequences) {
        for (Sequence seq : sequences) {
            if (seq.getTaxon().equals(id)) return seq;
        }
        return null;
    }

    public String toString() {
        return getTaxon() + ":" + getData();
    }


} // class Sequence
