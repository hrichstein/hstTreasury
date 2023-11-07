import os
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt

def mkPlot(x1,y1,x2,y2,lab1,lab2,outname='None',saveDir='./'):

    fig,ax = plt.subplots(figsize=(6,6))

    ax.scatter(x1,y1,label=f'{lab1}',s=20,color='black')
    ax.scatter(x2,y2,label=f'{lab2}',s=5,color='magenta')
    ax.set_rasterization_zorder(1)
    
    plt.savefig(os.path.join(saveDir,f'{outname}.pdf'),bbox_inches='tight',dpi=150)
    plt.close()

    return None


def matchIn(prime,cat,matchtol,xpri,ypri,xcat,ycat,idcat):

    """
    Input:
    prime: an ascii table, n by n, with the relevant information in columns
    cat: the array that will be searched for matches, n by n with relevant
    info in columns
    xpri: the column name where the x-coordinates are listed in the prime array
    ypri: the column name where the y-coordinates are listed in the prime array
    xcat: the column name where the x-coordinates are listed in the "cat"
    (matching) array
    ycat: the column name where the y-coordinates are listed in the "cat"
    (matching) array
    idcat: column in the matching array with the indices of sources as listed
    in the raw file

    Output:
    prime: a shortened version of the prime array that was input;
    only sources in the prime than have a match in the cat are listed
    match_ids: an len(prime) by 1 array that has the indices of the best-
    matching sources from the input cat array (ids come from the id column of
    the cat array; if the cat array was shortened from a longer cat array (to
    which this index list will be applied, the indices should have come from
    that longer cat array)). I'm matching based on index position, although in
    the future, to make things easier (possibly?) use real (kind of random,
    not necessarily listed order) IDs and match based on what the column ID
    says
    """

    matchids_in = np.zeros((len(prime),1))
    prime['dist'] = np.zeros((len(prime),1))

    nF = True
    row = 0

    while nF:
        # Finding the difference in the x and y positions between the prime
        # array, row by row
        # Checking that the difference in x and y (separately) is smaller than
        # the match tolerance
        # The statements in the brackets return indices for which the
        # conditions are true
        # Then, matchrows is a listing of the rows in the cat array where the
        # sources meet the positional requirements.
        matchrows = cat[(abs(prime[f'{xpri}'][row] - cat[f'{xcat}'])
                        <= matchtol) & (abs(prime[f'{ypri}'][row] - cat[f'{ycat}'])
                        <= matchtol)]

        # If only one source met the tolerance criteria, the index value for
        # that source is put into the matchids array. It will go in the row
        # corresponding to the row where the prime source was.
        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[f'{idcat}'][0]
            prime['dist'][row] = np.sqrt((prime[f'{xpri}'][row]
                                          - matchrows[0][f'{xcat}'])**2
                                          + (prime[f'{ypri}'][row]
                                          - matchrows[0][f'{ycat}'])**2)
            row += 1

        # If there is more than one source that meets the criteria,
        # I calculate the distance between all of the match sources
        # and the prime source. I put these in an array and find the (row)
        # index of the minimum distance. I proceed to put the cat array
        # source index into the matchids array.
        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt((prime[f'{xpri}'][row]
                                       - matchrows[f'{xcat}'][dd])**2
                                       + (prime[f'{ypri}'][row]
                                       - matchrows[f'{ycat}'][dd])**2)
            small = np.argsort(distDiff)
            idx = small[0]  # The index of the smallest distance in distDiff
            matchids_in[row][0] = matchrows[f'{idcat}'][idx]
            prime['dist'] = distDiff[idx] 
            row += 1

        # If there is nothing that meets the criteria, the prime source
        # row is removed, as well as the corresponding row in the matchids
        # array.
        else:
            prime.remove_row(row)
            matchids_in = np.delete(matchids_in,row,0)

        # If the row counter is longer than the length of the prime,
        # we've reached the end of the distance tabulations. I do a uniqueness
        # check, as if there's a repeat in the match_ids, it means multiple
        # prime sources matched with the same cat source. (I'm going
        # methodically through the prime list, but not removing best-matching
        # sources from the cat array.)
        if (row >= len(prime)):
            matchIDout, udx = np.unique(matchids_in,return_index=True)
            # udx is the array of unique indices. I check if this is less than
            # the prime length. If so, I use udx to get the relevant, unique
            # sources.
            if len(udx)<len(prime):
                uu = 1
                nF = False
                matchout = matchIDout
                # sorted unique matchID values
                # gets rid of duplicates if a source from "cat" matched multiple "prime" sources

            elif len(udx)==len(prime):
                uu = 0
                nF = False
                matchout = matchids_in
    return prime, matchout, uu
    

def matchlistID(prime,cat,matchtol,xpri,ypri,xcat,ycat,idcat,tag='',saveDir='./',plot=False):

    p_out, pids, uu = matchIn(prime,cat,matchtol,xpri,ypri,xcat,ycat,idcat)
    pidx = np.asarray(pids,dtype=int)  # Index values to CAT
    # p_out is the trimmed prime list, which has a column of indices that relate back to the original prime
    # pids are the index values of the cat array that correspond to the cat values that match up to something in p_out
    # if uu==0, it should be a one-to-one; if not, p_out would be longer than pids
    if plot:
        mkPlot(p_out[f'{xpri}'],p_out[f'{ypri}'],
               cat[pidx][f'{xcat}'],cat[pidx][f'{ycat}'],
               lab1='Prime1',lab2='Cat1',outname=f'aperRun1{tag}',
               saveDir=saveDir)
    
    if uu==0:
        primeOut = p_out
        match_ids = pids
        
    else:
        # Adding an id column to p_out, because this will be cut down
        # We'll need the indices that have a match in shortened cat to pull it out
        p_out['id_u'] = np.arange(0,len(p_out),1,dtype=int)

        # Length of prime sources with match in catalog (one cat source could match two or more prime sources)
        print(f'Number of matching prime sources: {len(p_out)}')  
        print(f'Unique cat values: {len(pids)}')  # Length of unique cat values that had matches in prime
        print(f'Length of original cat: {len(cat)}')  # Length of original cat
        
        print('Some overlap in matching indices')
        print(f'Original p_out length: {len(p_out):d}')
        p2 = cat[pidx]  # Getting only the unique cat values that had a match in p_out (shortened prime); 
        # Should have 'id' column
        c2 = deepcopy(p_out)  # Copying p_out to be put in as new "cat"

        print(f'Length of prime2, from cat: {len(p2)}')
        print(f'Length of cat2, from prime: {len(c2)}')

        # Running again with prime and cat inverted
        # Each value in the "cat" only appears once
        nF = True
        iter = 0
        while nF:
            if iter%2==0:
                print(f'ITER {iter:d}')
                print(f'Input length prime2, from cat: {len(p2)}')
                c_out, cids, uu = matchIn(p2,c2,matchtol,xcat,ycat,xpri,ypri,'id_u') #p2 is from CAT, c2 is from PRIME
                cidx = np.asarray(cids,dtype=int)
                print(f'Output length c_out, from cat: {len(c_out)}')
                print(f'Output length cids, from prime: {len(cids)}')
                if plot:
                    mkPlot(c_out[f'{xcat}'],c_out[f'{ycat}'],c2[f'{xpri}'],c2[f'{ypri}'],
                           lab1=f'Cat{iter}',lab2=f'Prime{iter}',outname=f'iter{iter}{tag}',
                           saveDir=saveDir)
            # c_out is the cut-down p2 (copy of cat) (if sources in the cut-down prime don't match the cut-down cat anymore)
            # cids are the index values of c2 (copy of p_out) that match with sources in c_out
                
            else:  # iter is odd, uu==1
                print(f'ITER {iter:d}')
                print(f'Input length prime2, from prime: {len(p2)}')
                c_out, cids, uu = matchIn(p2,c2,matchtol,xpri,ypri,xcat,ycat,'id_u')  #p2 is from PRIME, c2 is from CAT
                cidx = np.asarray(cids,dtype=int)
            # c_out is cut down p2 (from PRIME)
            # cids are index values of c2 (apply to c2 to get array with 'id' values?)
                sCatIdx = np.asarray(c2[cidx][idcat],dtype=int)  # Should be actual 'id' values from cat
                
                print(f'Output length of c_out, from prime: {len(c_out)}')
                print(f'Output length of cids, from cat: {len(cids)}')
                if plot:
                    mkPlot(c_out[f'{xpri}'],c_out[f'{ypri}'],c2[cidx][f'{xcat}'],
                           c2[cidx][f'{ycat}'],lab1=f'Prime{iter}',lab2=f'Cat{iter}',
                           outname=f'iter{iter}{tag}',saveDir=saveDir)

            if uu==0 and iter%2==0:  # uu has gotten to be 0 from the above if/else, even loop; can output and end
                print(f'ITER {iter:d}')
                primeOut = c2[cidx]  # c2 from prime here
                match_ids= c_out[f'{idcat}']  # Should be the original cat idx values that will match the primeOut list
                mcidx = np.asarray(match_ids,dtype=int)
                if plot:
                    mkPlot(primeOut[f'{xpri}'],primeOut[f'{ypri}'],cat[mcidx][f'{xcat}'],
                           cat[mcidx][f'{ycat}'],lab1=f'Prime{iter}',lab2=f'Cat{iter}',
                           outname=f'end{iter}{tag}',saveDir=saveDir)
                nF = False

            elif uu==1 and iter%2==0:  # If the problem hasn't been solved after the last iteration, which was even,
                # Send to the else loop immediately after while (switch prime, cat values below)
                print(f'ITER {iter:d}')
                p2 = deepcopy(c2[cidx])  # Individual prime values, 2 cuts deep, last copy of c2, which would be from PRIME
                c_out['id_u'] = np.arange(0,len(c_out),1,dtype=int)  # To be used in matching algorithm with c2 input
                c2 = deepcopy(c_out)  # Copy of cat, 2 cuts deep, with a two cat values that match to one prime value
                iter += 1

            elif uu==1:  # If still problem, but odd iteration, send to even loop, but switch prime and cat
                sCat = cat[sCatIdx]  # Shortened cat
                p2 = deepcopy(sCat)
                c2 = c_out
                c2['id_u'] = np.arange(0,len(c_out),1,dtype=int)
                iter += 1

            elif uu==0:  # If uu==0 but iter is odd; 
                primeOut = c_out
                match_ids = sCatIdx
                mcidx = np.asarray(match_ids,dtype=int)
                if plot:
                    mkPlot(primeOut[f'{xpri}'],primeOut[f'{ypri}'],cat[mcidx][f'{xcat}'],
                           cat[mcidx][f'{ycat}'],lab1=f'Prime{iter}',lab2=f'Cat{iter}',
                           outname=f'end{iter}{tag}',saveDir=saveDir)
                nF = False

                print(f'Length prime Out: {len(c_out)}')
                print(f'Length Cat IDs Out: {len(match_ids)}')

    return primeOut,match_ids