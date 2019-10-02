sequencePosteriors[sequence_, sitePrior_, siteProbs_, backgroundProbs_] :=
	Module[{smallerSequence=List[], sequencePosterior=List[], iteration=1, motifSize=Length[siteProbs]},
		Do[
			Do[AppendTo[smallerSequence, sequence[[j]]], {j, iteration, iteration+motifSize-1}];
			AppendTo[sequencePosterior, sitePosterior[smallerSequence, sitePrior, 1.0-sitePrior, siteProbs, backgroundProbs]];
			iteration=iteration+1;
			smallerSequence={},
		{i, Length[sequence]-motifSize+1}];
	sequencePosterior]
	
updateMotifPrior[normalizedPosteriors_] := 
	Module[{},
	    Sum[normalizedPosteriors[[i]][[j]], {i, Length[normalizedPosteriors]}, {j, Length[normalizedPosteriors[[i]]]}]/Sum[Length[normalizedPosteriors[[z]]],{z, Length[normalizedPosteriors]}]
	]
	
updatePFMCounts[motifLength_, input_, normalizedPosteriors_, motifPseudocounts_, erasers_] :=
	Module[{allNucleotides=List[], nucleotidesAtPositionX=List[], nucleotideStart=1, motifPosX=List[],
		PosFreqA=0, PosFreqC=0, PosFreqG=0, PosFreqT=0, position=1, sequence=1, PFM=List[], w, x},
		(*Runs through all sequences in our input, and resets nucleotideStart=1 at each iteration, since nucleotideStart
		is where we start iterating for each sequence in our input/*)
		Do[
			(*Runs thorugh the motif length, and appends nucleotidesAtPositionX to allNucleotides*)
			Do[
				(*Appends all nucleotides that can be of a certain position to nucleotidesAtPositionX*)
				Do[AppendTo[nucleotidesAtPositionX, input[[i]][[j]]],
				{j, nucleotideStart, Length[input[[i]]]-motifLength+1+(nucleotideStart-1)}];
			AppendTo[allNucleotides, nucleotidesAtPositionX];
			nucleotidesAtPositionX={};
			nucleotideStart++,
			{motifLength}];
		nucleotideStart=1,
		{i, 1, Length[input]}];
	    (*At the end of the 3 Do loops described above, allNucleotides will be of size Length[input]*motifLength, where each
	    entry in the list will be all nucleotides observed for a certain position, for a certain sequence in our input. For
	    example, if we are given a motifLength of 3, and the following input: {{4,3,4,2,3},{4,4,3,2,1}}, then allNucleotides[[1]]
	    will be {4,3,4} (all possible nucleotides that can take on position 1 from our first sequence in our input), 
	    allNucleotides[[2]] will be {3,4,2} (all possible nucleotides that can take on position 1 from our first sequence 
	    in our input), allNucleotides[[4]] will be {4,4,3,4} (possible nucleotides that can take on position 1 from our second 
	    sequence in our input). Notice that for a motif of size n (a1, a2, ....., ai), and an input of size p, all possible
	    nucleotides at position ai will be at {allNucleotides[[i]], allNucleotides[[i+p/i]], allNucleotides[[i+2(p/i)]], ...}*)
	    
	    (*Creates our position frequency matrix, initially all being the pseudocounts for each position in the motif.*)
	    Do[AppendTo[PFM, motifPseudocounts],{motifLength}];
	    
	    (*If erasers is false, we then have to initialize all of it's entries to 1.*)
	    (*If[erasers===False, erasers=Table[Table[1, Length[input[[isize]]]], {isize, Length[input]}]];*)
	    
	    (*Iterates once for each motif (to be added to the PFM)*)
	    Do[
	    	(*Sets the position frequencies for A, C, G, T based on the PFM we set above.*)
	    	PosFreqA=PFM[[position]][[1]]; PosFreqC=PFM[[position]][[2]]; PosFreqG=PFM[[position]][[3]]; PosFreqT=PFM[[position]][[4]];
	    	(*Runs through our allNucleotides list, in steps of motifLength (since we are trying to find all probabilities at a given position.)*)
	    	For[w=position, w<=Length[allNucleotides],
	    		(*Runs through all the nucleotides at that given position, in that specific sequence.*)
	    		For[x=1, x<=Length[allNucleotides[[w]]],
	    			(*Checks to see if erasers is false, if so we do not factor in the erasers.*)
	    			(*The if statements check to see if a nucleotide for a given position (w) in our sequence of possibilities (x) is
	    			A, C, G, T; it then adjusts by finding the normalized posterior relative to that sequence.*)
	    			If[erasers===False, 
	    				If[allNucleotides[[w]][[x]]==1, PosFreqA=PosFreqA+normalizedPosteriors[[sequence]][[x]],
	    					If[allNucleotides[[w]][[x]]==2, PosFreqC=PosFreqC+normalizedPosteriors[[sequence]][[x]],
	    						If[allNucleotides[[w]][[x]]==3, PosFreqG=PosFreqG+normalizedPosteriors[[sequence]][[x]], PosFreqT=PosFreqT+normalizedPosteriors[[sequence]][[x]]]]],
	    				If[allNucleotides[[w]][[x]]==1, PosFreqA=PosFreqA+normalizedPosteriors[[sequence]][[x]]*erasers[[sequence]][[x+position-1]],
	    					If[allNucleotides[[w]][[x]]==2, PosFreqC=PosFreqC+normalizedPosteriors[[sequence]][[x]]*erasers[[sequence]][[x+position-1]],
	    						If[allNucleotides[[w]][[x]]==3, PosFreqG=PosFreqG+normalizedPosteriors[[sequence]][[x]]*erasers[[sequence]][[x+position-1]], PosFreqT=PosFreqT+normalizedPosteriors[[sequence]][[x]]*erasers[[sequence]][[x+position-1]]]]]];
	    			x++;
	    			]; 
	    			w=w+motifLength; 
	    			sequence++;
	    		];
	    (*Sets the frequencies for a motif at position X*)
	    motifPosX={PosFreqA, PosFreqC, PosFreqG, PosFreqT};
	    (*Sets our probability frequency matrix at that position (which we have been keeping track thorugh the main outer loop, to
	    the recently calculated frequencies.*)
	    PFM[[position]]=motifPosX;
	    (*Resets frequencies to be calculated for next position in our motif.*)
	    motifPosX={};
	    position++;
	    (*Resets sequence to start at 1, to be able to find normalizedPosteriors.*)
	    sequence=1,
	    {motifLength}];
	    PFM
	]


sitePosterior[sequence_, sitePrior_, backgroundPrior_, siteProbs_, backgroundProbs_] := 
	Module[{i,likelihood},
		likelihood = Product[siteProbs[[i,sequence[[i]]]],{i,1,Length[sequence]}]*sitePrior; (*calculate likelyhood*)
		likelihood/(likelihood+backgroundPrior*Product[backgroundProbs[[sequence[[i]]]],{i,1,Length[sequence]}]) (*Calculate posterior probability*)
	]