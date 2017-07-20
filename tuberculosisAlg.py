# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 16:49:00 2017

@author: annaparker
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 20:42:57 2017

@author: annaparker
"""

import numpy as np
import random as rand

def readFile():
    """
    Reads the collection of strings Dna into a list
    """
    input = open("tubData.txt", "r");
    inp = []
    #read in input to text file
    for line in input:
        #print(line)
        inp += [line.strip('\n')]
    return inp



def SymbolToNumber(char): #turns given number into one of the letters of dna
    """Takes in a character found in DNA and turns it into a number"""
    if char == 'A':
        return 0
    elif char == 'C':
        return 1
    elif char == 'G':
        return 2
    elif char == 'T':
        return 3




def ProfileMostProbable(Text, k, Profile):
    """
    Given a string Text, finds the most probable kmer of length k given the profile matrix Profile
    
    """
    mostProbScore = 0
    mostProb = Text[0:k]
    
    #calculates the probability score for all possible kmers found in Text
    for i in range(0, len(Text)-k+1): 
        kmer = Text[i: i+k]
        currProbScore = 1
        
        #finds probability score by multiplying probabilities for each separate char based on profile matrix
        for j in range(0, k):
            currNum = SymbolToNumber(kmer[j])
            currProbScore = currProbScore * float(Profile[currNum][j])
                      
        #replaces mostProbScore if a higher probability score is found
        if currProbScore > mostProbScore:
            mostProbScore = currProbScore
            mostProb = kmer

    #returns kmer with highest probability             
    return mostProb  


  
def formProfileIgnore(K, kmers, ignore):
    """
    Makes a probability profile based on the kmers given of K length, not including the kmer in indexes ignore1 & ignore2.
    """    
    
    profile = np.ones((4, K)) #THIS STEP IMPLEMENTS PSEUDOCOUNTS 
    for i in range(0, len(kmers)): #do loop for each kmer except ignore1 & ignore2
        if(i not in ignore):
            currKmer = kmers[i] 
            for j in range(0, K): #for each character in kmer: increment that char in profile  
                currIndex = SymbolToNumber(currKmer[j]) 
                profile[currIndex][j] += 1
        

    #for each index in profile, divide value by the total number of kmers given
    for i in range(0, 4):
        for j in range(0, K):
            profile[i][j] = profile[i][j] / len(kmers)            
    
    
    return profile   





def Score(kmers, k, profile):
    """
    Calculates the similarity score of a group of strings kmer
    """
    
    #mostLikely is an array of the most probable character at each position in the kmer
    mostLikely = np.argmax(profile, axis = 0)
        
    
    totalScore = 0
    
    for i in range(0, k): #do loop for each position in kmer
        count = 0
        for j in range(0, len(kmers)): #do loop for number of strings in kmer
            currKmer = kmers[j]
            currIndex = SymbolToNumber(currKmer[i])
            if(currIndex != mostLikely[i]): #increment count each time a char doesn't match the most likely char
                count += 1
        
        totalScore += count #total score is aggregation of count for each position in kmer  
    return totalScore


def randomlySelect(k, t, dna):
    """
    creates a random list of motifs of length k from the t strings in dna
    """    
    probs = np.ones(len(dna[0]) - k + 1)
    probs = probs / np.sum(probs)

    motifList = []
    
    #for each string, randomly choose the first index from the probs buckets & add to motifList
    for i in range(0, t):
        firstIndex = np.digitize(np.random.random(), np.cumsum(probs))
        motif = dna[0][firstIndex: firstIndex+k]
        motifList += [motif]

    return motifList


def formProfile(K, kmers):
    """
    Makes a probability profile based on the kmers given of K length.
    """    
    
    profile = np.ones((4, K)) #THIS STEP IMPLEMENTS PSEUDOCOUNTS 
   
    for i in range(0, len(kmers)): #do loop for each kmer
        currKmer = kmers[i] 

        for j in range(0, K): #for each character in kmer: increment that char in profile  
            currIndex = SymbolToNumber(currKmer[j]) 
            profile[currIndex][j] += 1
    
    

    #for each index in profile, divide value by the total number of kmers given
    for i in range(0, 4):
        for j in range(0, K):
            profile[i][j] = profile[i][j] / len(kmers)            
    
    
    return profile 

def random(num):
    """
    returns a random number between 0 and num-1
    """
    probs = np.ones(num)
    probs = probs / np.sum(probs)
    x = np.digitize(np.random.random(), np.cumsum(probs))
    #print(x)    
    return x



def singleRandomProbableKmer(dna, k, profile):
    """
    based on profile returns the most probable kmer from within dna
    """
    #indices of array probs are the initial index of all possible kmers in dna
    probs = np.ones(len(dna) - k + 1) 

    #for each possible kmer in dna, calculate probability of each value in kmer & multiply together to put kmer probability in array prob 
    for i in range (0, len(dna) - k + 1):
        kmer = dna[i:i+k]        
        kmerProb = 1
        for j in range (0, k):
            val = SymbolToNumber(kmer[j])
            kmerProb *= profile[val][j]
            
        probs[i] = kmerProb 
    
    #normalize probs
    probs = probs / np.sum(probs)
    
    
    #choose kmer randomly from buckets sized depending on probability of kmer
    firstIndex = np.digitize(np.random.random(), np.cumsum(probs))
    kmer = dna[firstIndex: firstIndex+k]
    
    return kmer



def RandomizedMotifSearch(k, t, dna):
    """
    Given a collection of strings dna, t (the number of strings in dna) and k (the size of the k-mers), returns a collection of strings BestMotifs  
    """
    
    motifs = randomlySelect(k, t, dna) #randomly select starting motifs
    bestMotifs = motifs      
    profile = formProfile(k, bestMotifs)
    bestScore = Score(bestMotifs, k, profile)    
    
    
    forever = True
    while(forever): #run loop until it doesn't return a new best score
        
        profile = formProfile(k, motifs)        
        motifs = [] 
        
        for j in range(0, t): #for each dna string, find most likely kmer considering profile
            motifs += [ProfileMostProbable(dna[j], k, profile)] #update motifs with new kmers
        
        #update profile & calculate score given motif
        profile = formProfile(k, motifs)
        motifScore = Score(motifs, k, profile)        
        
        #if newly found score is greater than the best score, return the best score         
        if(motifScore >= bestScore):
            return bestScore, bestMotifs
        
        #otherwise, update bests & continue in the loop
        bestMotifs = motifs
        bestScore = motifScore



def GibbsSampler(dna, k, t, N, num_changes, startMotif):
    """
    given a collection of t strings dna, searches for the most similar kmers from each string. 
    
    Still works more cautiously than randomized motif search, but can be made`more aggressive than traditional Gibbs because it discards two
    kmers (instead of just one) from a current set of motifs at each iteration. Loop runs N times.
    """
   
    #randomly select the first set of motifs from dna
    motifs = startMotif   
    bestMotifs = startMotif
    profile = formProfile(k, bestMotifs)
    bestScore = Score(bestMotifs, k, profile)
    
    #in each loop, randomly select a dna string from which to select a new kmer probabalistically  
    for j in range(1, N):
        randStrings = rand.sample(range(t), num_changes)

        #make profile from motifs, not including the two motifs[randStrings]
        profile = formProfileIgnore(k, motifs, randStrings)
        
        for m in range(0, num_changes):                                
            motifs[randStrings[m]] = singleRandomProbableKmer(dna[randStrings[m]], k, profile)

        #update profile to calculate score        
        profile = formProfile(k, motifs)
        motifScore = Score(motifs, k, profile)
 

       #update bestScore if a better smotifScore is a better (smaller) score         
        if(motifScore < bestScore):
            bestMotifs = motifs
            bestScore = motifScore
            
    return bestScore, bestMotifs


def ShiftMotifs (BestRandomized, dna, k, numShift):
    """
    Takes in the collection of best motifs BestRandomized and shifts them over
    the number of times in numShift. F
    
    or example, if numshift is 1, returns the 
    original collection of BestRandomized, all strings in that collection shifted
    once, and all strings in that colletion shifted twice.  (The list of best motifs
    has tripled.)
    """
    shifted_motifs = []
    list_values = list(BestRandomized.values())
    shifted_motifs = list_values
    
    #shifting 1
    for j in range(1, numShift + 1):
        for motifs in list_values:
            current_set = []
            for i in range(len(motifs)):
                current_dna = dna[i]
                starting_index = current_dna.index(motifs[i])
                if((starting_index +k) < len(current_dna)):
                    current_set.append(current_dna[starting_index: starting_index + k])
                else: 
                    current_set.append(current_dna[starting_index + j: starting_index + k])
            shifted_motifs.append(current_set)    
    print(shifted_motifs)
    return shifted_motifs
        
 


                
def write(topScore, topMotif, randomizedBest, N, numChanges, numBestRandom): #write graph out to file
    target = open('tuberAlg.txt', 'w')
    target.write("N: %s\n" % N)
    target.write("numChanges %s\n" % numChanges)
    target.write("numBestRandom %s\n" % numBestRandom)
    target.write("topScore: %s\n" % topScore)
    target.write("topMotif:\n")
    for item in topMotif:
        target.write("%s\n" % item)
    target.write("\n \n \n")
    target.write("randomizedBest:\n")
    
    for score in randomizedBest:
        target.write("Score: %s\n" % score)
        for motifs in randomizedBest[score]: 
            target.write("%s\n" % motifs)
        target.write("\n")
    target.close()           





#####START OF MAIN

#input variables   
K = 20
T = 10
N = 100
numOfBestRandom = 50 #number of best motifs to collect from randomized


#POSSIBLE ALGORITHMIC CHANGES: 
#change gibbs aggression level:
numChanges = 5 #number of dna sequences to change each iteration of Gibbs
#shift best randomized motif:
numShift = 3 #number of times to shift motifs over


DNA = readFile()
topScore = np.inf

RandomizedMotifSearch(K, T, DNA)

best_randomized = {}

#RUN randomized motif search 1000 times, solution is best solutio of all runs
for i in range(0, 1000):
    best = RandomizedMotifSearch(K, T, DNA)
    
    #add the first ten motifs returned from randomized to the dictionary(keys are score values are motifs)
    if len(best_randomized) < numOfBestRandom:
        best_randomized[best[0]] = best[1]
    
    #if there are motifs returned from randomized that have a better score then return add those motifs to top_ten
    else: 
        max_score = max(list(best_randomized.keys()))#find max score in the dictionary
        
        if(best[0] < max_score): #if new best score found, update topScore & topMotif
            del best_randomized[max_score]
            best_randomized[best[0]] = best[1]

#RUN gibbs sampler once for each motif in the dictionary of bestMotifs found by RandomizedMotifeSearch
#topScore/ topMotifs hold the overall top values
for topMotif in best_randomized:
    best = GibbsSampler(DNA, K, T, N, numChanges, best_randomized[topMotif])

    if(best[0] < topScore): #if new best score found, update topScore & topMotif
        topScore = best[0]
        topMotifs = best[1]


write(topScore, topMotifs, best_randomized, N, numChanges, numOfBestRandom)
