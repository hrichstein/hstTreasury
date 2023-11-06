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


def matchIn(master,cat,matchtol,xmas,ymas,xcat,ycat,idcat):

    """
    Input:
    master: an ascii table, n by n, with the relevant information in columns
    cat: the array that will be searched for matches, n by n with relevant
    info in columns
    xmas: the column name where the x-coordinates are listed in the master array
    ymas: the column name where the y-coordinates are listed in the master array
    xcat: the column name where the x-coordinates are listed in the "cat"
    (matching) array
    ycat: the column name where the y-coordinates are listed in the "cat"
    (matching) array
    idcat: column in the matching array with the indices of sources as listed
    in the raw file

    Output:
    master: a shortened version of the master array that was input;
    only sources in the master than have a match in the cat are listed
    match_ids: an len(master) by 1 array that has the indices of the best-
    matching sources from the input cat array (ids come from the id column of
    the cat array; if the cat array was shortened from a longer cat array (to
    which this index list will be applied, the indices should have come from
    that longer cat array)). I'm matching based on index position, although in
    the future, to make things easier (possibly?) use real (kind of random,
    not necessarily listed order) IDs and match based on what the column ID
    says
    """

    matchids_in = np.zeros((len(master),1))
    master['dist'] = np.zeros((len(master),1))

    nF = True
    row = 0

    while nF:
        # Finding the difference in the x and y positions between the master
        # array, row by row
        # Checking that the difference in x and y (separately) is smaller than
        # the match tolerance
        # The statements in the brackets return indices for which the
        # conditions are true
        # Then, matchrows is a listing of the rows in the cat array where the
        # sources meet the positional requirements.
        matchrows = cat[(abs(master[f'{xmas}'][row] - cat[f'{xcat}'])
                        <= matchtol) & (abs(master[f'{ymas}'][row] - cat[f'{ycat}'])
                        <= matchtol)]

        # If only one source met the tolerance criteria, the index value for
        # that source is put into the matchids array. It will go in the row
        # corresponding to the row where the master source was.
        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[f'{idcat}'][0]
            master['dist'][row] = np.sqrt((master[f'{xmas}'][row]
                                          - matchrows[0][f'{xcat}'])**2
                                          + (master[f'{ymas}'][row]
                                          - matchrows[0][f'{ycat}'])**2)
            row += 1

        # If there is more than one source that meets the criteria,
        # I calculate the distance between all of the match sources
        # and the master source. I put these in an array and find the (row)
        # index of the minimum distance. I proceed to put the cat array
        # source index into the matchids array.
        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt((master[f'{xmas}'][row]
                                       - matchrows[f'{xcat}'][dd])**2
                                       + (master[f'{ymas}'][row]
                                       - matchrows[f'{ycat}'][dd])**2)
            small = np.argsort(distDiff)
            idx = small[0]  # The index of the smallest distance in distDiff
            matchids_in[row][0] = matchrows[f'{idcat}'][idx]
            master['dist'] = distDiff[idx] 
            row += 1

        # If there is nothing that meets the criteria, the master source
        # row is removed, as well as the corresponding row in the matchids
        # array.
        else:
            master.remove_row(row)
            matchids_in = np.delete(matchids_in,row,0)

        # If the row counter is longer than the length of the master,
        # we've reached the end of the distance tabulations. I do a uniqueness
        # check, as if there's a repeat in the match_ids, it means multiple
        # master sources matched with the same cat source. (I'm going
        # methodically through the master list, but not removing best-matching
        # sources from the cat array.)
        if (row >= len(master)):
            matchIDout, udx = np.unique(matchids_in,return_index=True)
            # udx is the array of unique indices. I check if this is less than
            # the master length. If so, I use udx to get the relevant, unique
            # sources.
            if len(udx)<len(master):
                uu = 1
                nF = False
                matchout = matchIDout
                # sorted unique matchID values
                # gets rid of duplicates if a source from "cat" matched multiple "master" sources

            elif len(udx)==len(master):
                uu = 0
                nF = False
                matchout = matchids_in
    return master, matchout, uu
    

def matchlistID(master,cat,matchtol,xmas,ymas,xcat,ycat,idcat,saveDir='./'):

    m_out, mids, uu = matchIn(master,cat,matchtol,xmas,ymas,xcat,ycat,idcat)
    midx = np.asarray(mids,dtype=int)  # Index values to CAT
    # m_out is the trimmed master list, which has a column of indices that relate back to the original master
    # mids are the index values of the cat array that correspond to the cat values that match up to something in m_out
    # if uu==0, it should be a one-to-one; if not, m_out would be longer than mids

    mkPlot(m_out[f'{xmas}'],m_out[f'{ymas}'],cat[midx][f'{xcat}'],cat[midx][f'{ycat}'],
           lab1='Master1',lab2='Cat1',outname='aperRun1',saveDir=saveDir)
    
    if uu==0:
        masterOut = m_out
        match_ids = mids
        
    else:
        # Adding an id column to m_out, because this will be cut down
        # We'll need the indices that have a match in shortened cat to pull it out
        m_out['id_u'] = np.arange(0,len(m_out),1,dtype=int)

        # Length of master sources with match in catalog (one cat source could match two or more master sources)
        print(f'Number of matching master sources: {len(m_out)}')  
        print(f'Unique cat values: {len(mids)}')  # Length of unique cat values that had matches in master
        print(f'Length of original cat: {len(cat)}')  # Length of original cat
        
        print('Some overlap in matching indices')
        print(f'Original m_out length: {len(m_out):d}')
        m2 = cat[midx]  # Getting only the unique cat values that had a match in m_out (shortened master); 
        # Should have 'id' column
        c2 = deepcopy(m_out)  # Copying m_out to be put in as new "cat"

        print(f'Length of master2, from cat: {len(m2)}')
        print(f'Length of cat2, from master: {len(c2)}')

        # Running again with master and cat inverted
        # Each value in the "cat" only appears once
        nF = True
        iter = 0
        while nF:
            if iter%2==0:
                print(f'ITER {iter:d}')
                print(f'Input length master2, from cat: {len(m2)}')
                c_out, cids, uu = matchIn(m2,c2,matchtol,xcat,ycat,xmas,ymas,'id_u') #m2 is from CAT, c2 is from MAS
                cidx = np.asarray(cids,dtype=int)
                print(f'Output length c_out, from cat: {len(c_out)}')
                print(f'Output length cids, from master: {len(cids)}')
                mkPlot(c_out[f'{xcat}'],c_out[f'{ycat}'],c2[f'{xmas}'],c2[f'{ymas}'],
                       lab1=f'Cat{iter}',lab2=f'Mas{iter}',outname=f'iter{iter}',saveDir=saveDir)
            # c_out is the cut-down m2 (copy of cat) (if sources in the cut-down master don't match the cut-down cat anymore)
            # cids are the index values of c2 (copy of m_out) that match with sources in c_out
                
            else:  # iter is odd, uu==1
                print(f'ITER {iter:d}')
                print(f'Input length master2, from master: {len(m2)}')
                c_out, cids, uu = matchIn(m2,c2,matchtol,xmas,ymas,xcat,ycat,'id_u')  #m2 is from MAS, c2 is from CAT
                cidx = np.asarray(cids,dtype=int)
            # c_out is cut down m2 (from MAS)
            # cids are index values of c2 (apply to c2 to get array with 'id' values?)
                sCatIdx = np.asarray(c2[cidx][idcat],dtype=int)  # Should be actual 'id' values from cat
                
                print(f'Output length of c_out, from master: {len(c_out)}')
                print(f'Output length of cids, from cat: {len(cids)}')
                mkPlot(c_out[f'{xmas}'],c_out[f'{ymas}'],c2[cidx][f'{xcat}'],c2[cidx][f'{ycat}'],
                       lab1=f'Mas{iter}',lab2=f'Cat{iter}',outname=f'iter{iter}',saveDir=saveDir)

            if uu==0 and iter%2==0:  # uu has gotten to be 0 from the above if/else, even loop; can output and end
                print(f'ITER {iter:d}')
                masterOut = c2[cidx]  # c2 from master here
                match_ids= c_out[f'{idcat}']  # Should be the original cat idx values that will match the masterOut list
                mcidx = np.asarray(match_ids,dtype=int)
                mkPlot(masterOut[f'{xmas}'],masterOut[f'{ymas}'],cat[mcidx][f'{xcat}'],cat[mcidx][f'{ycat}'],
                       lab1=f'Mas{iter}',lab2=f'Cat{iter}',outname=f'end{iter}',saveDir=saveDir)
                nF = False

            elif uu==1 and iter%2==0:  # If the problem hasn't been solved after the last iteration, which was even,
                # Send to the else loop immediately after while (switch master, cat values below)
                print(f'ITER {iter:d}')
                m2 = deepcopy(c2[cidx])  # Individual master values, 2 cuts deep, last copy of c2, which would be from MAS
                c_out['id_u'] = np.arange(0,len(c_out),1,dtype=int)  # To be used in matching algorithm with c2 input
                c2 = deepcopy(c_out)  # Copy of cat, 2 cuts deep, with a two cat values that match to one master value
                iter += 1

            elif uu==1:  # If still problem, but odd iteration, send to even loop, but switch master and cat
                sCat = cat[sCatIdx]  # Shortened cat
                m2 = deepcopy(sCat)
                c2 = c_out
                c2['id_u'] = np.arange(0,len(c_out),1,dtype=int)
                iter += 1

            elif uu==0:  # If uu==0 but iter is odd; 
                masterOut = c_out
                match_ids = sCatIdx
                mcidx = np.asarray(match_ids,dtype=int)
                mkPlot(masterOut[f'{xmas}'],masterOut[f'{ymas}'],cat[mcidx][f'{xcat}'],cat[mcidx][f'{ycat}'],
                       lab1=f'Mas{iter}',lab2=f'Cat{iter}',outname=f'end{iter}',saveDir=saveDir)
                nF = False

                print(f'Length Master Out: {len(c_out)}')
                print(f'Length Cat IDs Out: {len(match_ids)}')

    return masterOut,match_ids